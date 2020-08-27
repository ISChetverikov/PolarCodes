#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/GaussianApproximation.h"
#include "../include/ScFlipProgDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0


ScFlipProgDecoder::ScFlipProgDecoder(PolarCode * codePtr) : BaseDecoder(codePtr) {
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
	_crcPtr = new CRC(_codePtr->CrcPoly());
	_subchannelsMeansGa = std::vector<double>(n, 0);;
}

domain ScFlipProgDecoder::GetDomain() {
	return LLR;

}

void ScFlipProgDecoder::SetSigma(double sigma) {
	BaseDecoder::SetSigma(sigma);

	GaussianApproximation ga(sigma);
	size_t n = _codePtr->N();
	for (size_t i = 0; i < n; i++)
	{
		_subchannelsMeansGa[i] = ga.GetMu(i + 1, n);
	}

	return;
}

double ScFlipProgDecoder::f(double llr1, double llr2) {
	double sign = 1.0;

	if (llr1 < 0) {
		sign *= -1;
		llr1 *= -1;
	}
	if (llr2 < 0) {
		sign *= -1;
		llr2 *= -1;
	}

	return ((llr1 < llr2) ? llr1 : llr2) * sign;
}

double ScFlipProgDecoder::g(double llr1, double llr2, int u1) {
	return llr2 + (1 - 2 * u1) * llr1;
}

int ScFlipProgDecoder::L(double llr) {
	return (llr >= 0) ? 0 : 1;
}

void ScFlipProgDecoder::FillLeftMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<double>::iterator outIt,
	size_t n)
{
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++)
	{
		*outIt = f(*leftIt, *rightIt);
	}
}

void ScFlipProgDecoder::FillRightMessageInTree(std::vector<double>::iterator leftIt,
	std::vector<double>::iterator rightIt,
	std::vector<int>::iterator uhatIt,
	std::vector<double>::iterator outIt,
	size_t n)
{
	for (size_t i = 0; i < n; i++, leftIt++, rightIt++, outIt++, uhatIt++)
	{
		*outIt = g(*leftIt, *rightIt, *uhatIt);
	}
}


size_t ScFlipProgDecoder::log2(size_t n) {
	size_t m = 0;
	while (n > 0)
	{
		n = n >> 1;
		m++;
	}

	return m;
}

void ScFlipProgDecoder::PassDown(size_t iter) {
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

void ScFlipProgDecoder::PassUp(size_t iter) {
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

void ScFlipProgDecoder::DecodeFrom(int position) {
	size_t n = _codePtr->N();
	size_t m = _codePtr->m();

	PassUp(position);
	for (size_t i = position + 1; i < n; i++)
	{
		PassDown(i);
		if (_mask[i]) {
			_x[i] = L(_beliefTree[m][i]);
		}
		else {
			_x[i] = FROZEN_VALUE;
		}
		_uhatTree[m][i] = _x[i];
		PassUp(i);
	}
}

bool ScFlipProgDecoder::IsCrcPassed(std::vector<int> codeword) {
	size_t n = codeword.size();
	size_t k = _codePtr->k();
	size_t deg = _codePtr->CrcDeg();

	auto wordBits = _codePtr->UnfrozenBits();
	std::vector<int> word(n, 0);
	for (size_t i = 0; i < k; i++)
	{
		word[i] = codeword[wordBits[i]];
	}

	auto crcBits = _codePtr->CrcUnfrozenBits();
	std::vector<int> crc(deg, 0);
	for (size_t i = 0; i < deg; i++)
	{
		crc[i] = codeword[crcBits[i]];
	}

	auto crcReal = _crcPtr->Calculate(word);
	return crc == crcReal;
}

// mask: 1 - black, 0 - white
std::vector<int> ScFlipProgDecoder::GetCriticalSet(std::vector<int> mask, int position) {
	std::vector<std::vector<int>> tree;
	size_t m = _codePtr->m();
	size_t n = _codePtr->N();
	std::vector<int> result;

	for (size_t i = 0; i < m + 1; i++)
	{
		std::vector<int> temp((int)1 << i, 0);
		tree.push_back(temp);
	}
	for (size_t i = 0; i < n; i++)
	{
		tree[m][i] = mask[i];
	}

	for (int i = (int)m - 1; i >= 0 ; i--)
	{
		size_t levelSize = tree[i].size();
		for (int j = 0; j < levelSize; j++)
		{
			int left = j << 1;
			int right = left + 1;

			tree[i][j] = tree[i + 1][left] && tree[i + 1][right];
		}
	}

	for (size_t i = 0; i < m + 1; i++)
	{
		int levelSize = (int)tree[i].size();
		for (int j = 0; j < levelSize; j++)
		{
			if (!tree[i][j])
				continue;

			int criticalBit = j << (m - i);
			if (criticalBit > position)
				result.push_back(criticalBit);

			int left = j;
			int right = j;
			for (size_t k = i + 1; i <= m; i++)
			{
				left = left << 1;
				right = (right << 1) + 1;

				for (size_t l = left; l <= right; l++)
					tree[k][l] = 0;
			}
		}
	}

	return result;

}

// with using means after gaussian approxiamtion procedure
std::vector<int> ScFlipProgDecoder::SortCriticalBits(std::vector<int> criticalSet, std::vector<double> llrs) {
	size_t length = criticalSet.size();
	std::vector<int> result(length, 0);
	std::vector<double> sortingMetric(length, 0);
	
	for (size_t i = 0; i < length; i++)
	{
		int criticalInd = criticalSet[i];
		sortingMetric[i] = fabs(llrs[criticalInd]) / _subchannelsMeansGa[criticalInd];
	}

	for (size_t i = 0; i < length; i++)
	{
		auto minIt = std::min_element(sortingMetric.begin(), sortingMetric.end());
		int minInd = (int)std::distance(sortingMetric.begin(), minIt);
		result[i] = criticalSet[minInd];

		criticalSet.erase(criticalSet.begin() + minInd);
		sortingMetric.erase(sortingMetric.begin() + minInd);
	}

	return result;
}

std::vector<int> ScFlipProgDecoder::Decode(std::vector<double> inLlr) {
	size_t n = inLlr.size();
	size_t m = _codePtr->m();
	size_t k = _codePtr->k();
	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = inLlr[i];
	}

	for (size_t i = 0; i < n; i++)
	{

		PassDown(i);
		if (_mask[i]) {
			_x[i] = L(_beliefTree[m][i]);
		}
		else {
			_x[i] = FROZEN_VALUE;
		}
		_uhatTree[m][i] = _x[i];
		PassUp(i);
	}

	if (!IsCrcPassed(_x)) {
		auto criticalSet = GetCriticalSet(_mask, 0);
		std::vector<int> suspectedBits = SortCriticalBits(criticalSet, _beliefTree[m]);
		for (size_t i = 0; i < suspectedBits.size(); i++)
		{
			int bitPosition = suspectedBits[i];
			_x[bitPosition] = !_x[bitPosition];

			DecodeFrom(bitPosition);

			if (IsCrcPassed(_x))
				break;
		}
	}

	std::vector<int> result(k, 0);
	size_t l = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (_mask[i] == 0)
			continue;

		result[l] = _x[i];
		l++;
	}

	return result;
}