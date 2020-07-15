#include <vector>
#include <map>
#include <unordered_map>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScFanoDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"
#include "../include/GaussianApproximation.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScFanoDecoder::ScFanoDecoder(PolarCode * codePtr, double T, double delta) : BaseDecoder(codePtr) {
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();
	_treeHeight = m + 1;

	for (size_t i = 0; i < _treeHeight; i++)
	{
		std::vector<double> b(n, -10000.0);
		std::vector<int> u(n, -1);

		_beliefTree.push_back(b);
		_uhatTree.push_back(u);
	}

	_mask = _codePtr->BitsMask();
	_x = std::vector<int>(n, -1);
	_T = T;
	_delta = delta;
}

domain ScFanoDecoder::GetDomain() {
	return P1;
}


double ScFanoDecoder::f(double p1, double p2) {

	return p1 * (1 - p2) + p2 * (1 - p1);
}

double ScFanoDecoder::g(double p1, double p2) {
	if (p1 == 0 && p2 == 1 || p2 == 0 && p1 == 1)
		return 1 / 2;

	return p1 * p2 / (p1 * p2 + (1 - p1) * (1 - p2));
}

int ScFanoDecoder::L(double p1) {
	return (p1 >= 1 / 2) ? 1 : 0;
}

void ScFanoDecoder::FillLeftMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<double>::iterator outIt,
	size_t n) 
{
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++)
	{
		*outIt = f(*leftIt, *rightIt);
	}
}

void ScFanoDecoder::FillRightMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<int>::iterator uhatIt,
	std::vector<double>::iterator outIt,
	size_t n)
{
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++, uhatIt++)
	{
		auto r = f(*uhatIt, *leftIt);
		*outIt = g(r, *rightIt);
	}
}


size_t log2(size_t n) {
	size_t m = 0;
	while (n > 0)
	{
		n = n >> 1;
		m++;
	}

	return m;
}

void ScFanoDecoder::PassDown(size_t iter) {
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();

	size_t iterXor;
	size_t level;
	if (iter) {
		iterXor = iter ^ (iter - 1);
		level = m - log2(iterXor);
	}
	else {
		level = 0;
	}

	std::vector<int> binaryIter(m - level, 0);
	size_t iterCopy = iter;
	for (int i = (int)binaryIter.size() - 1; i >= 0; i--)
	{
		binaryIter[i] = iterCopy % 2;
		iterCopy = iterCopy >> 1;
	}
	
	size_t length = (size_t)1 << (m - level - 1);
	for (size_t i = level; i < m; i++)
	{
		size_t ones = ~0u;
		size_t offset = iter & (ones << (m - i));
		
		if (!binaryIter[i - level]) {
			FillLeftMessageInTree(_beliefTree[i].begin() + offset,
				_beliefTree[i].begin() + offset + length,
				_beliefTree[i + 1].begin() + offset,
				length);
		}
		else {
			
			FillRightMessageInTree(_beliefTree[i].begin() + offset,
				_beliefTree[i].begin() + offset + length,
				_uhatTree[i + 1].begin() + offset,
				_beliefTree[i + 1].begin() + offset + length,
				length);
		}
		
		length = length / 2;
	}
}

void ScFanoDecoder::PassUp(size_t iter) {
	size_t m = _codePtr->m();
	size_t iterCopy = iter;

	size_t bit = iterCopy % 2;
	size_t length = 1;
	size_t level = m;
	size_t offset = iter;
	while (bit != 0)
	{
		offset -= length;
		for (size_t i = 0; i < length; i++)
		{
			_uhatTree[level - 1][offset + i] = _uhatTree[level][offset + i] ^ _uhatTree[level][offset + length + i];
		}
		for (size_t i = 0; i < length; i++)
		{
			_uhatTree[level - 1][offset + length + i] = _uhatTree[level][offset + length + i];
		}
		

		iterCopy = iterCopy >> 1;
		bit = iterCopy % 2;
		length *= 2;
		level -= 1;
	}
}

void UpdateT(double & T, double & delta, double & tau) {
	while (T + delta < tau)
		T += delta;
}

void BackwardMove(std::vector<double> & beta, std::vector<bool> & gamma, std::vector<int> & mask, int & i, double & T, double & delta, bool & B, size_t & firstInfoBit) {

	while (true) {

		double mu = 0;
		if (i != 0) {
			mu = beta[i - 1];
		}

		if (mu >= T && i != firstInfoBit) {
			bool gammaPrevious = gamma[i];
			// Find the previous info bit
			i--;
			while (!mask[i])
				i--;

			if (!gamma[i]) {
				B = true;
				return;
			}
		}
		else {
			B = false;
			T = T - delta;
			return;
		}
	}

	return;

}

std::vector<int> ScFanoDecoder::Decode(std::vector<double> inP1) {
	size_t n = inP1.size();
	size_t m = _codePtr->m();
	size_t k = _codePtr->k();
	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inP1[i];
	}

	GaussianApproximation ga(_sigma);
	std::vector<double> p(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}

	int i = 0;
	bool B = false;
	std::vector<double> beta(n, 0.0);
	std::vector<bool> gamma(n, 0);
	double T = _T;
	double delta = _delta;

	size_t firstInfoBit = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (_mask[i]) {
			firstInfoBit = i;
			break;
		}
	}
	while (i < n)
	{
		PassDown(i); // get p1 metric in _beliefTree[m][i]
		if (_mask[i]) {

			double p0 = 1 - _beliefTree[m][i];
			double previous = (i == 0) ? 0 : beta[i - 1];
			double m0 = previous + log(p0 / (1 - p[i]));

			double p1 = _beliefTree[m][i];
			double m1 = previous + log(p1 / (1 - p[i]));

			double max = (m1 > m0) ? m1 : m0 ;
			int argmax = (m1 > m0) ? 1 : 0;

			double min = (m1 > m0) ? m0 : m1;
			int argmin = (m1 > m0) ? 0 : 1;

			if (max > T) {
				if (!B) {
					_x[i] = argmax;
					_uhatTree[m][i] = _x[i];
					PassUp(i);
					beta[i] = max;
					gamma[i] = false;

					double mu = 0;
					if (i != 0) {
					//if (i != firstInfoBit) {
						/*int previousInfoBit = i;
						previousInfoBit--;
						while (!_mask[previousInfoBit])
							previousInfoBit--;*/
						mu = beta[i-1];
					}

					if (mu < T + delta) 
						UpdateT(T, delta, mu);
					i++;
				}
				else {
					if (min > T) {
						_x[i] = argmin;
						_uhatTree[m][i] = _x[i];
						PassUp(i);
						beta[i] = min;
						gamma[i] = true;
						B = false;
						i++;
					}
					else {
						if (i == firstInfoBit) {
							T = T - delta;
							B = false;
						}
						else {
							BackwardMove(beta, gamma, _mask, i, T, delta, B, firstInfoBit);
						}
					}
				}
			}
			else {
				if (i == firstInfoBit) {
					T = T - delta;
				}
				else
					BackwardMove(beta, gamma, _mask, i, T, delta, B, firstInfoBit);
			}
		}
		else {
			_x[i] = 0;
			_uhatTree[m][i] = _x[i];
			PassUp(i);
			if (i != 0) {
				beta[i] = beta[i - 1];
				gamma[i] = gamma[i - 1];
			}
			else {
				beta[i] = 0;
				gamma[i] = 0;
			}
				
			
			i++;
		}
	}
	
	// Get info bits
	std::vector<int> result(k, 0);
	size_t j = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (_mask[i] == 0)
			continue;

		result[j] = _x[i];
		j++;
	}

	return result;
}