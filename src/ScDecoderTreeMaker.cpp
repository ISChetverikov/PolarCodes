#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/ScDecoderTreeMaker.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"
#include "../include/GaussianApproximation.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScDecoderTreeMaker::ScDecoderTreeMaker(PolarCode * codePtr, double approximationSigma) : ScDecoder(codePtr) {
	size_t n = _codePtr->N();
	GaussianApproximation ga(approximationSigma);
	_p = std::vector<double>(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		_p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}
}

std::vector<int> ScDecoderTreeMaker::Decode(std::vector<double> inP1) {
	_path = "";

	size_t n = inP1.size();
	size_t m = _codePtr->m();
	size_t k = _codePtr->k();
	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inP1[i];
	}

	int i = 0;
	int j = -1;
	bool B = false;
	std::vector<double> beta(k, 0.0); // only for frozen bits
	std::vector<double> metrics(n, 0); // for all bits
	std::vector<bool> gamma(k, 0);
	std::vector<int> A = _codePtr->UnfrozenBits(); // info set
	
	size_t firstInfoBit = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (_maskWithCrc[i]) {
			firstInfoBit = i;
			break;
		}
	}
	while (i < n)
	{
		PassDown(i); // get p1 metric in _beliefTree[m][i]
		double p0 = 1 - _beliefTree[m][i];
		double p1 = _beliefTree[m][i];

		if (_maskWithCrc[i]) {
			double previous = (i == 0) ? 0 : metrics[i - 1];
			double m0 = previous + log(p0 / (1 - _p[i]));
			double m1 = previous + log(p1 / (1 - _p[i]));

		}
		else {
			_x[i] = FROZEN_VALUE;
			_uhatTree[m][i] = _x[i];
			PassUp(i);

			double currentMetric = log(p0) - log(1 - _p[i]);
			// cumulative
			metrics[i] = currentMetric;
			if (i != 0)
				metrics[i] += metrics[i - 1];

			i++;
		}

		// HERE trace
		// _path += std::to_string(i) + ": " + std::to_string(metrics[i]) + "\n";
		///////
	}

	return TakeResult();
}