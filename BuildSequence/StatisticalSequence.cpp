#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <mpi.h>
#include <ctime>
#include <chrono>
#include <sstream>
#include <iomanip>

#include "../include/ScDecoder.h"
#include "../include/PolarCode.h"
#include "../include/Encoder.h"
#include "../include/BpskBscChannel.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/Exceptions.h"

#include "StatisticalSequence.h"

using std::vector;
using std::chrono::steady_clock;
using std::chrono::time_point;
using std::cout;

vector<vector<int>> GenerateTrialFrozenSets(int n, vector<int> leader, int ProcRank, int ProcNum) {
	
	vector<int> remainedBits;

	for (size_t i = 0; i < n; i++)
	{
		if (std::find(leader.begin(), leader.end(), i) != leader.end())
			continue;

		remainedBits.push_back(i);
	}


	vector<vector<int>> result;
	for (size_t j = ProcRank; j < remainedBits.size(); j += ProcNum)
	{
		vector<int> copy = leader;
		copy.push_back(remainedBits[j]);
		result.push_back(copy);
	}

	int max_count = remainedBits.size() / ProcNum;
	if (remainedBits.size() % ProcNum != 0)
		max_count += 1;
	
	//cout << "Rank: " << ProcRank << " max_count: " << max_count << std::endl;

	if (result.size() < max_count)
		result.push_back(result[result.size() - 1]);

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
	double sigma_right = 100.0;
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

	//double snr = pow(2, 2 * capacity) - 1 - 10 * log10(R);
	double snr = -10 * log10(2 * sigma *sigma) - 10 * log10(R);

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

std::string GetDumpFilename(int k, int n, double snr) {
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	std::ostringstream oss;
	oss << std::setprecision(3);
	oss << snr;
	oss << "_";
	oss << std::put_time(&tm, "%m-%d_%H-%M-%S");
	auto snrTimeStr = oss.str();

	std::string dumpFilename =  "leader_" + std::to_string(n) +
		"_" + std::to_string(k) + "_" + snrTimeStr + ".dump";

	return dumpFilename;
}

std::string GetOuputFilename(int n, int k, double snr) {
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	std::ostringstream oss;
	oss << std::setprecision(3);
	oss << snr;
	oss << "_";
	oss << std::put_time(&tm, "%m-%d-%H-%M-%S");
	auto snrTimeStr = oss.str();

	std::string dumpFilename = "sequence_" + std::to_string(n) +
		"_" + std::to_string(k) + "_"  + snrTimeStr + ".dump";

	return dumpFilename;
}

// Only the root process writes into the files of the folder
void BuiltSequenceStatistically(std::string folder, int m, int k, int maxTestsCount, int maxRejectionsCount, double _snr,
	int ProcRank, int ProcNum) {

	size_t n = 1 << m;
	
	vector<int> leader = {};

	time_point<steady_clock> t1;
	std::string dumpFilename;
	if (ProcRank == 0) {

		std::cout << "SNR is omitted: " << _snr << std::endl;

		t1 = steady_clock::now();
		dumpFilename = folder + "leader.dump";//GetDumpFilename(k, n, snr);

		TryCreateDump(dumpFilename);
		TryLoadDump(dumpFilename, leader);
		//std::reverse(leader.begin(), leader.end());
		
	}
	int leaderSize = leader.size();
	MPI_Bcast(&leaderSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank != 0)
		leader = vector<int>(leaderSize);
	MPI_Bcast(leader.data(), leaderSize, MPI_INT, 0, MPI_COMM_WORLD);
	/*
	for (size_t i = 0; i < leader.size(); i++)
	{
		cout << leader[i] << " ";
	}
	cout << std::endl;*/

	// Init objects

	size_t k_current = leader.size() + 1;
	for (; k_current <= k; k_current++)
	{
		double snr = PickUpSnr(m, k_current);
		struct {
			double value;
			int   index;
		} p_best_struct;

		double p_best = 1.1;
		double p = 1.1;

		short break_flag = 0;
		vector<vector<int>> frozenCombinations;

		frozenCombinations = GenerateTrialFrozenSets(n, leader, ProcRank, ProcNum);
		
		/*int rank = 0;
		while (rank < ProcNum) {

			if (rank == ProcRank) {
				cout << "Rank: " << ProcRank << std::endl;
				for (size_t i = 0; i < frozenCombinations.size(); i++)
				{
					cout << "\t";
					for (size_t j = 0; j < frozenCombinations[i].size(); j++)
					{
						cout << frozenCombinations[i][j] << " ";
					}

					cout << std::endl;
				}

			}
			
			rank++;
			MPI_Barrier(MPI_COMM_WORLD);
		}*/
		
		for (int j = (int)frozenCombinations.size() - 1; j >= 0 ; j--)
		{

			PolarCode * codePtr = new PolarCode(m, frozenCombinations[j]);

			Encoder * encoderPtr = new Encoder(codePtr);
			BaseDecoder * decoderPtr = new ScDecoder(codePtr);
			BaseChannel * channelPtr = new BpskBscChannel();

			BaseSimulator * simulatorPtr = new MonteCarloSimulator(maxTestsCount, maxRejectionsCount, codePtr, encoderPtr, channelPtr, decoderPtr, 0);

			p = simulatorPtr->Run(snr).fer;
			//p_best_struct.index = ProcRank;
			
			//cout << "Rank: " << ProcRank << " p_value: " << p_best_struct.value << std::endl;
			
			//MPI_Allreduce(MPI_IN_PLACE, &p_best_struct, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

			//cout << "Rank: " << ProcRank << " p_best: " << p_best << std::endl;
			//cout << "Rank: " << ProcRank << " p_value: " << p_best_struct.value << " p_index: " << p_best_struct.index << std::endl;

			if (p < p_best) {
				
				leader = frozenCombinations[j];
				p_best = p;

				//int leaderSize = leader.size();
				//MPI_Bcast(&leaderSize, 1, MPI_INT, p_best_struct.index, MPI_COMM_WORLD);

				//MPI_Bcast(leader.data(), leader.size(), MPI_INT, p_best_struct.index, MPI_COMM_WORLD);
				//p_best = p_best_struct.value;
				//MPI_Bcast(&p_best, 1, MPI_DOUBLE, p_best_struct.index, MPI_COMM_WORLD);
				if (p_best == 0.0)
					break_flag = 1;
			}
			
			MPI_Allreduce(MPI_IN_PLACE, &break_flag, 1, MPI_SHORT, MPI_LOR, MPI_COMM_WORLD);

			delete simulatorPtr;
			delete channelPtr;
			delete decoderPtr;
			delete encoderPtr;
			delete codePtr;

			if (break_flag) {
				break;
			}
		}

		p_best_struct.value = p_best;
		p_best_struct.index = ProcRank;

		//cout << "Rank: " << ProcRank << " leaderSize: " << leader.size() << std::endl;

		MPI_Allreduce(MPI_IN_PLACE, &p_best_struct, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
		MPI_Bcast(leader.data(), leader.size(), MPI_INT, p_best_struct.index, MPI_COMM_WORLD);

		//MPI_Barrier(MPI_COMM_WORLD);
		if (ProcRank == 0) {
			std::cout << "------- " << k_current << " ------" << std::endl;
			std::cout << "SNR: " << snr << std::endl;
			std::cout << "P best: " << p_best << std::endl;
			std::cout << "-------------------------------" << std::endl;

			SaveDump(dumpFilename, leader);
		}

	}

	//if (ProcRank == 0) {
		//ClearDump(dumpFilename);
		//std::cout << "Dump is succesfully cleared..." << std::endl;
	//}
	
	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == 0) {
		auto t2 = std::chrono::steady_clock::now();
		auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		// continue barrier
		std::cout << "Time: " << dt << std::endl;

		std::string filename = folder + GetOuputFilename(n, k, _snr);
		std::ofstream file(filename);

		vector<int> sequence = bitsToSequence(n, leader);
		for (size_t i = 0; i < sequence.size(); i++)
		{
			file << sequence[i] << " ";
		}
	}
	

}