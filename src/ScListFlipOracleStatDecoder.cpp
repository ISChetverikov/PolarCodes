
#include <algorithm>
#include <sstream>
#include "../include/ScListFlipOracleStatDecoder.h"

#define FROZEN_VALUE 0

ScListFlipOracleStatDecoder::ScListFlipOracleStatDecoder(PolarCode * codePtr, int L) : ScListFlipStatDecoder(codePtr, L) {
}

int ScListFlipOracleStatDecoder::GetFirstErrorPositionMetric(std::vector<std::vector<int>> candidates, std::vector<double> metrics, std::vector<int> originalCodeword) {
	
	auto maxIt = std::max_element(metrics.begin(), metrics.end());
	int maxInd = (int)std::distance(metrics.begin(), maxIt);
	
	int i = 0;
	for (; i < originalCodeword.size(); i++)
	{
		if (candidates[maxInd][i] != originalCodeword[i])
			break;
	}

	return i;
}

int ScListFlipOracleStatDecoder::GetFirstErrorPositionMaxDistance(std::vector<std::vector<int>> candidates, std::vector<int> originalCodeword) {

	int n = (int)originalCodeword.size();

	int best = 0;
	for (int i = 0; i < _L; i++)
	{
		int j = 0;
		for (; j < n; j++)
		{
			if (candidates[i][j] != originalCodeword[j])
				break;
		}
		
		if (j > best)
			best = j;
	}

	return best;
}

void ScListFlipOracleStatDecoder::ClearStatistic() {
	_singleOracleFlipsStat.clear();
	_doubleOracleFlipsStat.clear();
	_tripleOracleFlipsStat.clear();

	_countErrorneous = 0;
	_singleSuccessfulFlips = 0;
	_doubleSuccessfulFlips = 0;
	_tripleSuccessfulFlips = 0;
}

std::string ScListFlipOracleStatDecoder::GetStatistic() {
	std::stringstream ss;
	
	/*ss << "Single Flip:\n";
	for (std::pair<int, int> single : _singleOracleFlipsStat)
	{
		ss << "(" << single.first << "): " << single.second << "\n";
	}

	ss << "Double Flip:\n";
	for (std::pair<std::tuple<int, int>, int> double_ : _doubleOracleFlipsStat)
	{
		ss << "(" << std::get<0>(double_.first) << ", " << std::get<1>(double_.first) << "): " << double_.second << "\n";
	}

	ss << "Triple Flip:\n";
	for (std::pair<std::tuple<int, int, int>, int> triple : _tripleOracleFlipsStat)
	{
		ss << "(" << std::get<0>(triple.first) << ", " << std::get<1>(triple.first) << ", " << std::get<2>(triple.first) << "): " << triple.second << "\n";
	}*/

	ss << "Count:" << _count << "\n";
	ss << "List Error Count: " << _countErrorneous << "\n";
	ss << "Single: " << _singleSuccessfulFlips << "\n";
	ss << "Double: " << _doubleSuccessfulFlips << "\n";
	ss << "Triple: " << _tripleSuccessfulFlips << "\n";
	
	return ss.str();
}

std::vector<int>  ScListFlipOracleStatDecoder::Decode(std::vector<double> beliefs) {
	std::vector<int> actualCodeword(beliefs.size(), 0);
	std::vector<int> flipPositions;
	bool isError = false;
	
	_count++;

	DecodeListInternal(beliefs);
	std::vector<double> metrics = _metrics; // Take result changes the metrics
	std::vector<int> result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError)
		return result;

	_countErrorneous++;

	// Single
	int firstErrorInd = GetFirstErrorPositionMaxDistance(_candidates, _codeword);
	flipPositions.clear();
	flipPositions.push_back(firstErrorInd);

	DecodeFlipListInternal(beliefs, flipPositions);

	metrics = _metrics;

	result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError) {
		_singleSuccessfulFlips++;
		if (_singleOracleFlipsStat.find(firstErrorInd) == _singleOracleFlipsStat.end())
			_singleOracleFlipsStat[firstErrorInd] = 1;
		else
			_singleOracleFlipsStat[firstErrorInd]++;

		return result;
	}
	

	// Doulbe
	int secondErrorInd = GetFirstErrorPositionMaxDistance(_candidates, _codeword);
	flipPositions.push_back(secondErrorInd);

	DecodeFlipListInternal(beliefs, flipPositions);

	metrics = _metrics;

	result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError) {
		_doubleSuccessfulFlips++;
		auto tuple = std::make_tuple(firstErrorInd, secondErrorInd);

		if (_doubleOracleFlipsStat.find(tuple) == _doubleOracleFlipsStat.end())
			_doubleOracleFlipsStat[tuple] = 1;
		else
			_doubleOracleFlipsStat[tuple]++;

		return result;
	}
	
	// Triple
	int thirdErrorInd = GetFirstErrorPositionMaxDistance(_candidates, _codeword);
	flipPositions.push_back(thirdErrorInd);

	DecodeFlipListInternal(beliefs, flipPositions);
	result = TakeListStatResult(actualCodeword);
	isError = actualCodeword != _codeword;

	if (!isError) {
		_tripleSuccessfulFlips++;
		auto tuple = std::make_tuple(firstErrorInd, secondErrorInd, thirdErrorInd);

		if (_tripleOracleFlipsStat.find(tuple) == _tripleOracleFlipsStat.end())
			_tripleOracleFlipsStat[tuple] = 1;
		else
			_tripleOracleFlipsStat[tuple]++;

		return result;
	}

	return result;
}