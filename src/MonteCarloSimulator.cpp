#include <random>
#include <string>
#include <fstream>
#include <chrono>
#include <exception>
#include "../include/OperationsCount.h"
#include "../include/SimulationParameters.h"
#include "../include/MonteCarloSimulator.h"
#include "../include/ScFanoDecoder.h"
#include "../include/ScDecoder.h"
#include "../include/ScListDecoder.h"
#include "../include/ScStackDecoder.h"

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

template<class T>
std::string VecToStr(std::vector<T> vec) {
	std::ostringstream streamObj;
	std::string result;

	for (size_t i = 0; i < vec.size(); i++)
	{
		streamObj << vec[i] << ", ";
		//result += std::to_string(vec[i]) + ", ";
	}

	return streamObj.str();
}

#ifdef PARALLEL_DECODER
std::vector<int> ReadSequenceFromFileParallel(std::string path) {
	std::vector<int> seq;
	std::string line;
	std::ifstream myFile(path);

	std::getline(myFile, line);

	int val;
	std::stringstream ss(line);

	while (ss >> val)
		seq.push_back(val);

	return seq;
}
#endif // PARALLEL_DECODER
SimulationIterationResults MonteCarloSimulator::Run(double snr)
{	
	SimulationIterationResults result;

	size_t n = _codePtr->N();
	size_t k = _codePtr->k();
	std::vector<int> word(k, 0);
	std::vector<int> codeword(n, 0);
	std::vector<double> channelOuput(n, 0);
	std::vector<int> decoded(n, 0);
	
	auto t1 = std::chrono::steady_clock::now();

	double sigma = GetSigma(snr, (double)k / n);
	int tests = 0;
	int wrong_dec = 0;
	//int wrong_bits = 0;

	std::random_device randomDevice;
	std::uniform_int_distribution<> uniform_discrete_dist(0, 1);

#ifdef PARALLEL_DECODER
	// HERE parallel decoder
	std::vector<int> parallelDecoded;
	size_t m = _codePtr->m();
	PolarCode * parallelCodePtr = new PolarCode(m, k, ReadSequenceFromFileParallel("C:\\Users\\ische\\source\\repos\\PolarCodes\\polar_sequences\\" + std::to_string(n) + ".txt"), _codePtr->CrcPoly());
	ScListDecoder parallelDecoder(parallelCodePtr, 4);
	//////
#endif // PARALLEL_DECODER

	_decoderPtr->SetSigma(sigma);
	_channelPtr->SetSnr(snr);

	_decoderPtr->ClearOperationsCount();

	while ((tests < _maxTestsCount || _maxTestsCount == -1) && (wrong_dec < _maxRejectionsCount)) {
		tests++;

		std::generate(word.begin(), word.end(), [&]() { return uniform_discrete_dist(randomDevice); });

		codeword = _encoderPtr->Encode(word);
		// Give answer to a decoder for debugging or statistic retreving
		_decoderPtr->SetDecoderAnswer(_encoderPtr->PolarTransform(codeword));

		channelOuput = _channelPtr->Pass(codeword);
		//channelOuput = { 3.2307, -3.11039, -3.36801, -6.56322, 2.06576, 1.72037, 2.46457, 2.86076, -0.109397, -1.24382, -4.95555, 4.86218, 1.07354, 2.6011, -2.87971, -3.12015 };

		decoded = _decoderPtr->Decode(channelOuput);

#ifdef PARALLEL_DECODER1
		// HERE parallel decoder to comparision, SC - worse decoder
		parallelDecoded = parallelDecoder.Decode(channelOuput);
		
		if (word != decoded && parallelDecoded == word) {
			std::cout << "Find" << std::endl;
			std::string debugInfo = _decoderPtr->GetPathInfo();
			std::string filename = "C:\\Users\\ische\\source\\repos\\PolarCodes\\results\\dump_SCS.debug";
			DumpInfo(filename, VecToStr<double>(channelOuput));
			DumpInfo(filename, VecToStr<int>(word));
			DumpInfo(filename, VecToStr<int>(decoded));
			DumpInfo(filename, VecToStr<int>(parallelDecoded));
			DumpInfo(filename, debugInfo);
			DumpInfo(filename, "");
		}
		///////////
#endif //PARALLEL_DECODER

		if (tests % 10000 == 0)
			std::cout << tests << std::endl;
		
		if (decoded != word) {
			wrong_dec += 1;

			/*for (size_t i = 0; i < n; i++)
			{
				if (decoded[i] != word[i])
					wrong_bits += 1;
			}*/
		}
	}

	//std::cout << "BER: " << (double)wrong_bits / n / tests << std::endl;

	OperationsCount operationsCount = _decoderPtr->GetOperationsCount();

	auto t2 = std::chrono::steady_clock::now();

	result.snr = snr;
	result.ebn0 = GetEbN0(snr, k, n);
	result.sigma = sigma;

	result.fer = (double)wrong_dec / tests;

	result.elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	result.testsCount = tests;
	result.rejectionsCount = wrong_dec;
	
	result.operationsCount = operationsCount;

	return result;
}
