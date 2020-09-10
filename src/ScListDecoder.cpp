#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../include/ScListDecoder.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScListDecoder::ScListDecoder(PolarCode * codePtr, int L) : ScCrcAidedDecoder(codePtr) {
	_L = L;

	size_t m = _codePtr->m();
	size_t n = _codePtr->N();
	_treeHeight = m + 1;

	for (size_t i = 0; i < _L; i++)
	{
		std::vector<std::vector<double>> beliefTree;
		std::vector<std::vector<int>> uhatTree;

		std::vector<double> b(n, -10000.0);
		std::vector<int> u(n, -1);

		for (size_t i = 0; i < _treeHeight; i++)
		{
			beliefTree.push_back(b);
			uhatTree.push_back(u);
		}

		_beliefTrees.push_back(beliefTree);
		_uhatTrees.push_back(uhatTree);
	}

	std::vector<int> x(n, -1);
	for (size_t i = 0; i < _L; i++)
	{
		_candidates.push_back(x);
	}
	_metrics = std::vector<double>(_L, 0);
	
	// optimization allocation
	_binaryIter = std::vector<int>(m, 0);
}

double ScListDecoder::StepMetric(double belief, int decision) {
#ifdef DOMAIN_LLR
	return (decision) ? -belief : belief;
#elif DOMAIN_P1
	return log((decision ^ (belief > 0.5))  ? 1 - belief : belief);
#endif // DOMAIN
}

void ScListDecoder::PassDownList(size_t iter) {
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

	//std::vector<int> binaryIter(m - level, 0);
	int size = (int)(m - level);
	size_t iterCopy = iter;
	for (int i = size - 1; i >= 0; i--)
	{
		_binaryIter[i] = iterCopy % 2;
		iterCopy = iterCopy >> 1;
	}

	size_t length = (size_t)1 << (m - level - 1);
	for (size_t i = level; i < m; i++)
	{
		size_t ones = ~0u;
		size_t offset = iter & (ones << (m - i));

		for (size_t j = 0; j < 2 * _L; j++)
		{
			if (!_binaryIter[i - level]) {
				FillLeftMessageInTree(_beliefTrees[j][i].begin() + offset,
					_beliefTrees[j][i].begin() + offset + length,
					_beliefTrees[j][i + 1].begin() + offset,
					length);
			}
			else {

				FillRightMessageInTree(_beliefTrees[j][i].begin() + offset,
					_beliefTrees[j][i].begin() + offset + length,
					_uhatTrees[j][i + 1].begin() + offset,
					_beliefTrees[j][i + 1].begin() + offset + length,
					length);
			}

		}
		
		length = length / 2;
	}
}

void ScListDecoder::PassUpList(size_t iter) {
	size_t m = _codePtr->m();
	size_t iterCopy = iter;

	size_t bit = iterCopy % 2;
	size_t length = 1;
	size_t level = m;
	size_t offset = iter;
	while (bit != 0)
	{
		offset -= length;
		for (size_t j = 0; j < 2 * _L; j++)
		{
			for (size_t i = 0; i < length; i++)
			{
				_uhatTrees[j][level - 1][offset + i] = _uhatTrees[j][level][offset + i] ^ _uhatTrees[j][level][offset + length + i];
			}
			for (size_t i = 0; i < length; i++)
			{
				_uhatTrees[j][level - 1][offset + length + i] = _uhatTrees[j][level][offset + length + i];
			}
		}
		
		iterCopy = iterCopy >> 1;
		bit = iterCopy % 2;
		length *= 2;
		level -= 1;
	}
}

void ScListDecoder::DecodeListInternal(std::vector<double> inLlr) {
	size_t n = inLlr.size();
	size_t m = _codePtr->m();

	// Fill each tree in the forrest with input llrs
	for (size_t j = 0; j < _L; j++)
		for (size_t i = 0; i < n; i++)
			_beliefTrees[j][0][i] = inLlr[i];
	
	int logL = log2(_L) - 1;
	// first log(_L) bits
	size_t i = 0; // number of bit (j - number of condidate in list)
	while (i < logL)
	{
		PassDownList(i);

		if (_maskWithCrc[i]) {
			int value = 0;
			for (size_t j = 0; j < _L; j++) {
				_candidates[j][i] = value;
				if (value % (i + 1) == 0)
					value = !value;
			}
		}
		else {
			for (size_t j = 0; j < _L; j++)
				_candidates[j][i] = FROZEN_VALUE;
		}

		for (size_t j = 0; j < _L; j++) {
			_uhatTrees[j][m][i] = _candidates[j][i];
			_metrics[j] += StepMetric(_beliefTrees[j][m][i], _candidates[j][i]);
		}
			
		PassUpList(i);

		if (_maskWithCrc[i])
			i++;
	}
	while (i < n)
	{
		PassDownList(i);

		if (_maskWithCrc[i]) {
			// 2 x 2 array
			auto sortedIndices = GetIndices();
			
			//for (size_t j = 0; j < _L; j++)
				//_candidates[j][i] = L(_beliefTrees[j][m][i]);
		}
		else {
			for (size_t j = 0; j < _L; j++)
				_candidates[j][i] = FROZEN_VALUE;
		}

		for (size_t j = 0; j < _L; j++)
			_uhatTrees[j][m][i] = _candidates[j][i];
		PassUpList(i);

		i++;
	}

	return;
}

std::vector<int> ScListDecoder::Decode(std::vector<double> inLlr) {

	DecodeInternal(inLlr);

	return TakeResult();
}