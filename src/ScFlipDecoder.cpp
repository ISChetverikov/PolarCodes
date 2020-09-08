#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/PolarCode.h"
#include "../include/ScFlipDecoder.h"
#include "../include/ScCrcAidedDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0


ScFlipDecoder::ScFlipDecoder(PolarCode * codePtr, domain domain, bool isMinSum, int T) : ScCrcAidedDecoder(codePtr, domain, isMinSum) {
	_T = T;
}

std::vector<int> ScFlipDecoder::GetSmallestLlrsIndices(std::vector<double> llrs, int count) {
	std::vector<int> indices(count, 0);
	for (size_t i = 0; i < llrs.size(); i++)
	{
		llrs[i] = fabs(llrs[i]);
	}

	for (size_t i = 0; i < count; i++)
	{
		auto minIt = std::min_element(llrs.begin(), llrs.end());
		auto minInd = (int)std::distance( llrs.begin(), minIt);
		indices[i] = minInd;
		llrs.erase(minIt);
	}
	
	return indices;
}

std::vector<int> ScFlipDecoder::Decode(std::vector<double> inLlr) {

	DecodeEnternal(inLlr);

	size_t m = _codePtr->m();
	if (_T > 0 && !IsCrcPassed(_x)) {
		std::vector<int> suspectedBits = GetSmallestLlrsIndices(_beliefTree[m], _T);
		for (size_t i = 0; i < _T; i++)
		{
			int bitPosition = suspectedBits[i];
			_x[bitPosition] = !_x[bitPosition];

			DecodeFrom(bitPosition);

			if (IsCrcPassed(_x))
				break;
		}
	}

	return TakeResult();
}