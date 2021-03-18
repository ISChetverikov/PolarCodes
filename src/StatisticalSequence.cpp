#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <omp.h>

#include "../include/StatisticalSequence.h"
#include "../include/ScDecoder.h"
#include "../include/PolarCode.h"
#include "../include/Encoder.h"
#include "../include/BpskBscChannel.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/Exceptions.h"

using std::vector;

vector<vector<int>> GenerateTrialFrozenSets(int n, vector<int> leader) {
	vector<vector<int>> result;

	for (size_t i = 0; i < n; i++)
	{
		if (std::find(leader.begin(), leader.end(), i) != leader.end())
			continue;

		vector<int> copy = leader;
		copy.push_back(i);
		result.push_back(copy);
	}

	return result;

}

vector<int> bitsToSequence(int n, vector<int> usedBits) {

	size_t k = usedBits.size();
	vector<int> result(n - k);

	for (size_t i = 0, j = 0; i < n; i++)
	{
		if (std::find(usedBits.begin(), usedBits.end(), i) != usedBits.end())
			continue;

		result[j] = i;
		j++;
	}


	result.insert(result.end(), usedBits.rbegin(), usedBits.rend());

	return result;

}

// mean = 1.0
double NormalCdfOfZero(double sigma) {
	return 1.0 / 2 * (1 + erf(-1.0 / sigma / sqrt(2)));
}

double BinaryEntropy(double p) {
	if (p == 0.0 || p == 1.0)
		return 1;

	return -p * log2(p) - (1 - p) * log2(1 - p);
}

double PickUpSnr(int m, int k) {

	int n = 1 << m;
	double R = (double) k / n;

	double eps = 1e-7;

	//double p = 0.5 - 0.5 * sqrt(2 - exp(x));
	
	// inverse capacity chanel dihotomy
	double p_left = 0;
	double p_right = 0.5;
	double p = 0.25;
	double capacity = 1 - BinaryEntropy(p);

	while ((p_right - p_left) >= eps) {
		p = p_left + (p_right - p_left) / 2;

		capacity = 1 - BinaryEntropy(p);

		if (capacity < R)
			p_right = p;
		else
			p_left = p;
	}

	// Normal CDF dihotomy relative to sigma
	double sigma_right = 2.0;
	double sigma_left = 0.1;
	double sigma = 1.0;
	double cdf_zero = NormalCdfOfZero(sigma);

	while ((sigma_right - sigma_left) >= eps) {
		sigma = sigma_left + (sigma_right - sigma_left) / 2;

		cdf_zero = NormalCdfOfZero(sigma);

		if (p < cdf_zero)
			sigma_right = sigma;
		else
			sigma_left = sigma;
	}

	double snr = -10 * log10(2 * sigma *sigma);

	return snr;
}

void TryCreateDump(std::string filename) {
	std::ofstream file(filename, std::ios::app);

	if (!file.is_open())
		throw FileIsNotOpennedException("Could not create dump file: " + filename);
}

void TryLoadDump(std::string filename, vector<int> & futureLeader) {
	std::ifstream file(filename);

	if (!file.is_open())
		throw FileIsNotOpennedException("Could not open for reading dump file: " + filename);

	int num;
	while (file >> num)
	{
		futureLeader.push_back(num);
	}
}

void SaveDump(std::string filename, vector<int> leader) {
	std::ofstream file(filename, std::ios::trunc);

	if (!file.is_open())
		throw FileIsNotOpennedException("Could not open for writing dump file: " + filename);

	for (size_t i = 0; i < leader.size(); i++)
	{
		file << leader[i] << " ";
	}
}

void ClearDump(std::string filename) {
	std::ofstream file(filename, std::ios::trunc);

	if (!file.is_open())
		throw FileIsNotOpennedException("Could not clear dump file: " + filename);
}

void BuiltSequenceStatistically(std::string folder, int m, int k, int maxTestsCount, int maxRejectionsCount) {

	size_t n = 1 << m;
	
	auto t1 = std::chrono::steady_clock::now();

	double snr =PickUpSnr(m, k);
	std::cout << "SNR: " << snr << std::endl;
	
	vector<int> leader = {};

	std::string dumpFilename = "leader.dump";

	TryCreateDump(dumpFilename);
	TryLoadDump(dumpFilename, leader);

	for (size_t k_current = leader.size() + 1; k_current <= k; k_current++)
	{
		vector<vector<int>> frozenCombinations = GenerateTrialFrozenSets(n, leader);

		double p_best = 1.0;

		volatile bool break_flag = false;

		#pragma omp parallel for shared(break_flag) num_threads(8)
		for (int j = (int)frozenCombinations.size() - 1; j >= 0 ; j--)
		{
			if (break_flag)
				continue;

			PolarCode * codePtr = new PolarCode(m, frozenCombinations[j]);

			Encoder * encoderPtr = new Encoder(codePtr);
			BaseDecoder * decoderPtr = new ScDecoder(codePtr);
			BaseChannel * channelPtr = new BpskBscChannel(double(k_current) / n);

			BaseSimulator * simulatorPtr = new MonteCarloSimulator(maxTestsCount, maxRejectionsCount, codePtr, encoderPtr, channelPtr, decoderPtr, 0);

			double p = simulatorPtr->Run(snr).fer;
			//std::cout << p << std::endl;
			
			#pragma omp critical
			{
				if (p <= p_best) {
					p_best = p;
					leader = frozenCombinations[j];
					if (p_best == 0.0)
						break_flag = true;
				}
			}

			delete simulatorPtr;
			delete decoderPtr;
			delete encoderPtr;
			delete codePtr;
		}

		std::cout << "P best: " << p_best << std::endl;
		std::cout << "------- " << k_current << " ------" << std::endl;

		SaveDump(dumpFilename, leader);
	}

	ClearDump(dumpFilename);
	std::cout << "Dump is succesfully cleared..." << std::endl;

	auto t2 = std::chrono::steady_clock::now();
	auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

	std::cout << "Time: " << dt << std::endl;

	std::string filename = folder + std::to_string(n);
	std::ofstream file(filename);
	//for (size_t i = 0; i < leader.size(); i++)
	//{
		//std::cout << leader[i] << " ";
		//file << leader[i] << " ";
	//}
	//std::cout << "\n";

	//file << std::endl;

	vector<int> sequence = bitsToSequence(n, leader);
	for (size_t i = 0; i < sequence.size(); i++)
	{
		file << sequence[i] << " ";
	}

}