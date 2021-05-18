#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm>

#include "../include/ScStackDecoder.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0
#define MIN_METRIC -100000.0

using std::pair;
using std::vector;

ScStackDecoder::ScStackDecoder(PolarCode * codePtr, int L, int D) : ScCrcAidedDecoder(codePtr) {
	_L = L;
	_D = D;

	size_t m = _codePtr->m();
	size_t n = _codePtr->N();
	_treeHeight = m + 1;

	std::vector<double> b(n, -10000.0);
	std::vector<int> u(n, -1);

	std::vector<std::vector<double>> beliefTree;
	std::vector<std::vector<int>> uhatTree;

	for (size_t i = 0; i < _treeHeight; i++)
	{
		beliefTree.push_back(b);
		uhatTree.push_back(u);
	}

	for (size_t i = 0; i < _D; i++)
	{
		_beliefTrees.push_back(beliefTree);
		_uhatTrees.push_back(uhatTree);
	}

	for (size_t i = 0; i < _D; i++)
	{
		_candidates.push_back(u);
	}

	_metrics = std::vector<pair<double, int>>(0);
	_current_bits = std::vector<int>(_D, 0);
	_paths_limits = std::vector<int>(n, 0);
	
}

double ScStackDecoder::StepMetric(double belief, int decision) {
#ifdef DOMAIN_LLR

#ifdef MINSUM

	return (belief < 0 && decision == 0 || belief > 0 && decision == 1) ? -fabs(belief) : 0;
#else
	const double limit = 10000000.0;

	if (belief < -limit)
		belief = -limit;

	if (belief > limit)
		belief = limit;

	double p0_methric = -log(1 + exp(-belief));
	double p1_methric = -log(1 + exp(belief));
	return (decision) ? p1_methric : p0_methric;

#endif // MINSUM

#elif DOMAIN_P1
	return log((decision) ? belief : 1 - belief);
#endif // DOMAIN
}

void ScStackDecoder::PassDownAll(size_t iter) {
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();

	size_t iterXor;
	size_t level;
	if (iter) {
		iterXor = iter ^ (iter - 1);
		level = m - FirstBitPos(iterXor);
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

		for (size_t j = 0; j < _D; j++)
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

void ScStackDecoder::PassDownStackLeader(size_t iter) {
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();

	size_t iterXor;
	size_t level;
	if (iter) {
		iterXor = iter ^ (iter - 1);
		level = m - FirstBitPos(iterXor);

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

		int j = _metrics[0].second;
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

		length = length / 2;

	}
}


void ScStackDecoder::PassUpAll(size_t iter) {
	size_t m = _codePtr->m();
	size_t iterCopy = iter;

	size_t bit = iterCopy % 2;
	size_t length = 1;
	size_t level = m;
	size_t offset = iter;
	while (bit != 0)
	{
		offset -= length;
		for (size_t j = 0; j < _D; j++)
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

void ScStackDecoder::PassUpStackSelectevely(size_t iter, vector<int> elements) {
	size_t m = _codePtr->m();
	size_t iterCopy = iter;

	size_t bit = iterCopy % 2;
	size_t length = 1;
	size_t level = m;
	size_t offset = iter;

	while (bit != 0)
	{
		offset -= length;
		for (const auto& j : elements)
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

void ScStackDecoder::EliminatePaths(size_t iter) {
	for (size_t i = 0; i < _D; i++)
	{
		int j = _metrics[i].second;
		if (_current_bits[j] <= iter) {
			_metrics[i].first = MIN_METRIC;
		}
	}

	std::sort(_metrics.rbegin(), _metrics.rend());
}


vector<int> ScStackDecoder::Decode(std::vector<double> inLlr) {
	size_t n = inLlr.size();
	size_t m = _codePtr->m();

	_metrics.clear();
	// Fill each tree in the forrest with input llrs
	for (size_t j = 0; j < _D; j++) {
		for (size_t i = 0; i < n; i++)
			_beliefTrees[j][0][i] = inLlr[i];

		_metrics.push_back(std::make_pair(0.0, j));
		_current_bits[j] = 0;
	}

	for (size_t i = 0; i < n; i++)
	{
		_paths_limits[i] = 0;
	}

	int logD = (int)FirstBitPos(_D) - 1;
	size_t i_all = 0; // number of bit (j - index of candidate in the stack)
	size_t i_unfrozen = 0;

	while (i_unfrozen < logD)
	{
		PassDownAll(i_all);

		if (_maskWithCrc[i_all]) {
			int value = 0;
			for (size_t j = 0; j < _D; j++) {
				_candidates[j][i_all] = value;
				_uhatTrees[j][m][i_all] = value;
				_metrics[j].first += StepMetric(_beliefTrees[j][m][i_all], _candidates[j][i_all]);
				_current_bits[j]++;

				if ((j + 1) % (1 << i_unfrozen) == 0) // all paths at the logD first steps 
					value = !value;
			}
			i_unfrozen++;
		}
		else {
			for (size_t j = 0; j < _D; j++) {
				_uhatTrees[j][m][i_all] = FROZEN_VALUE;
				_candidates[j][i_all] = FROZEN_VALUE;
				_current_bits[j] ++;
				_metrics[j].first += StepMetric(_beliefTrees[j][m][i_all], FROZEN_VALUE);
			}
		}

		PassUpAll(i_all);

		i_all++;
	}
	for (size_t i = 0; i < i_all - 1; i++)
	{
		_paths_limits[i] = _L;
	}
	
	std::sort(_metrics.rbegin(), _metrics.rend());

	bool isExited = false;
	while (!isExited)
	{
		if (_metrics[0].first == MIN_METRIC)
			break;

		i_all = _current_bits[_metrics[0].second];
		_paths_limits[i_all-1]++;
		if (_paths_limits[i_all-1] >= _L)
			EliminatePaths(i_all-1);

		PassDownStackLeader(i_all);
		std::vector<int> newElements(0);

		int leaderIndex = _metrics[0].second;
		int lastIndex = _metrics.back().second;

		if (_maskWithCrc[i_all]) {

			double step0 = StepMetric(_beliefTrees[leaderIndex][m][i_all], 0);
			double step1 = StepMetric(_beliefTrees[leaderIndex][m][i_all], 1);

			double metricsNew0 = _metrics[0].first + step0;
			double metricsNew1 = _metrics[0].first + step1;

			int highBit, lowBit;
			double highMetric, lowMetric;
			if (metricsNew1 > metricsNew0) {
				highBit = 1;
				lowBit = 0;
				highMetric = metricsNew1;
				lowMetric = metricsNew0;
			}
			else {
				highBit = 0;
				lowBit = 1;
				highMetric = metricsNew0;
				lowMetric = metricsNew1;
			}

			// high bit structures
			_candidates[leaderIndex][i_all] = highBit;
			_uhatTrees[leaderIndex][m][i_all] = highBit;
			_current_bits[leaderIndex]++;
			//_paths_limits[i_all]++;
		
			_metrics[0].first = highMetric;

			newElements.push_back(leaderIndex);
			// low bit structures
			if (lowMetric > _metrics.back().first) {
				_beliefTrees[lastIndex] = _beliefTrees[leaderIndex];
				_candidates[lastIndex] = _candidates[leaderIndex];
				_uhatTrees[lastIndex] = _uhatTrees[leaderIndex];
				_current_bits[lastIndex] = _current_bits[leaderIndex]; // already increased in highbit routine

				_candidates[lastIndex][i_all] = lowBit;
				_uhatTrees[lastIndex][m][i_all] = lowBit;
				_metrics.back().first = lowMetric;
				//_paths_limits[i_all]++;
				
				newElements.push_back(lastIndex);

			}

		}
		else {
			_candidates[leaderIndex][i_all] = FROZEN_VALUE;
			_uhatTrees[leaderIndex][m][i_all] = FROZEN_VALUE;
			newElements.push_back(leaderIndex); // leader remained
			_current_bits[leaderIndex]++;
			_metrics[0].first += StepMetric(_beliefTrees[leaderIndex][m][i_all], FROZEN_VALUE);

		}
		std::sort(_metrics.rbegin(), _metrics.rend());
		PassUpStackSelectevely(i_all, newElements);

		if (i_all == (n - 1)) {
			if (IsCrcPassed(_candidates[_metrics[0].second]) || _paths_limits[i_all] >= _L)
				isExited = true;
			else if (_paths_limits[i_all] < _L) {
				_paths_limits[i_all]++;
				_metrics[0].first = MIN_METRIC;
				std::sort(_metrics.rbegin(), _metrics.rend());
			}

		}
	}

	std::vector<int> result(_codePtr->k(), 0);
	std::vector<int> codewordBits = _codePtr->UnfrozenBits();
	for (size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = _candidates[_metrics[0].second][codewordBits[i]];
	}

	return result;
}

