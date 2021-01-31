#include <random>
#include <string>
#include <fstream>
#include <chrono>
#include <exception>
#include "../include/SimulationParameters.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/ScFanoDecoder.h"
#include "../include/ScDecoder.h"

MonteCarloSimulator::MonteCarloSimulator(int maxTestsCount,
	int maxRejectionsCount,
	PolarCode * codePtr,
	Encoder * encoderPtr,
	BaseChannel * channelPtr,
	BaseDecoder * decoderPtr,
	bool isSigmaDependOnR) : BaseSimulator(codePtr, encoderPtr, channelPtr, decoderPtr, isSigmaDependOnR)
{
	_maxTestsCount = maxTestsCount;
	_maxRejectionsCount = maxRejectionsCount;
}

void DumpInfo(std::string filename, std::string info) {
	std::ofstream resultsFileStream;
	resultsFileStream.open(filename, std::fstream::out | std::fstream::app);
	resultsFileStream << info << std::endl;

	resultsFileStream.close();
}

std::string VecToStr(std::vector<double> vec) {
	std::string result;

	for (size_t i = 0; i < vec.size(); i++)
	{
		result += std::to_string(vec[i]) + ", ";
	}

	return result;
}

SimulationIterationResults MonteCarloSimulator::Run(double snr)
{	
	SimulationIterationResults result;

	size_t n = _codePtr->N();
	size_t k = _codePtr->k();
	std::vector<int> word(k, 0);
	std::vector<int> codeword(n, 0);
	std::vector<double> channelOuput(n, 0);
	std::vector<int> decoded(n, 0);
	std::vector<int> parallelDecoded(n, 0);

	auto t1 = std::chrono::steady_clock::now();

	double sigma = GetSigma(snr, (double)k / n);
	int tests = 0;
	int wrong_dec = 0;
	
	std::random_device randomDevice;
	std::uniform_int_distribution<> uniform_discrete_dist(0, 1);

	ScDecoder parallelDecoder(_codePtr);

	_decoderPtr->SetSigma(sigma);
	_channelPtr->SetSnr(snr);

	while ((tests < _maxTestsCount || _maxTestsCount == -1) && (wrong_dec < _maxRejectionsCount)) {
		tests++;

		std::generate(word.begin(), word.end(), [&]() { return uniform_discrete_dist(randomDevice); });

		codeword = _encoderPtr->Encode(word);
		// Give answer to a decoder for debugging or statistic retreving
		_decoderPtr->SetDecoderAnswer(_encoderPtr->PolarTransform(codeword));

		channelOuput = _channelPtr->Pass(codeword);
		
		decoded = _decoderPtr->Decode(channelOuput);

		// HERE parallel decoder to comparision, SC - worse decoder
		parallelDecoded = parallelDecoder.Decode(channelOuput);

		if (word == decoded && parallelDecoded != word) {
			std::string debugInfo = _decoderPtr->GetPathInfo();
			std::string filename = "C:\\Users\\ische\\source\\repos\\PolarCodes\\results\\dump.debug";
			DumpInfo(filename, VecToStr(channelOuput));
			DumpInfo(filename, debugInfo);
			DumpInfo(filename, "");
		}
		///////////

		if (tests % 1000 == 0)
			std::cout << tests << std::endl;
		
		if (decoded != word)
			wrong_dec += 1;
	}

	auto t2 = std::chrono::steady_clock::now();

	result.snr = snr;
	result.ebn0 = GetEbN0(snr, k, n);
	result.sigma = sigma;

	result.fer = (double)wrong_dec / tests;

	result.elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	result.testsCount = tests;
	result.rejectionsCount = wrong_dec;
	
	return result;
}
