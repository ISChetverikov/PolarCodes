#include <sstream>
#include <algorithm>

#include "../include/ScCrcAidedDecoder.h"
#include "../include/ScFlipStatDecoder.h"

ScFlipStatDecoder::ScFlipStatDecoder(PolarCode * codePtr) : ScCrcAidedDecoder(codePtr) {
	_unfrozenBits = _codePtr->UnfrozenBits();
	ClearStatistic();
}

void ScFlipStatDecoder::ClearStatistic() {
	size_t k = _codePtr->k();
	_singleFlipStatistic = std::vector<int>(k, 0);
	_doubleFlipStatistic = std::vector<std::vector<int>>(k, _singleFlipStatistic);
}

std::string ScFlipStatDecoder::GetStatistic() {
	std::stringstream ss;
	std::sort(_singleFlipStatistic.rbegin(), _singleFlipStatistic.rend());
	ss << "Single Flip:\n";
	for (size_t i = 0; i < _codePtr->k(); i++)
	{
		if (_singleFlipStatistic[i])
			ss << "(" << _unfrozenBits[i] << "): " << _singleFlipStatistic[i] << "\n";
	}

	ss << "Double Flip:\n";
	for (size_t i = 0; i < _codePtr->k(); i++)
	{
		for (size_t j = 0; j < _codePtr->k(); j++)
		{
			std::sort(_doubleFlipStatistic[i].rbegin(), _doubleFlipStatistic[i].rend());
			if (_doubleFlipStatistic[i][j]) {
				ss << "(" << _unfrozenBits[i] << ", " << _unfrozenBits[j] << "): " << _doubleFlipStatistic[i][j] << "\n";
			}
		}
	}

	return ss.str();
}

std::vector<int>  ScFlipStatDecoder::Decode(std::vector<double> beliefs) {
	DecodeInternal(beliefs);
	
	auto originalCodeword = _codeword;

	if (_x == originalCodeword)
		return TakeResult();

	size_t k = _codePtr->k();
	for (size_t i = 0; i < k; i++)
	{
		int bitPosition = _unfrozenBits[i];

		if (_x[bitPosition] == originalCodeword[bitPosition])
			continue;

		_x[bitPosition] = !_x[bitPosition];

		DecodeFrom(bitPosition);

		if (_x == _codeword) {
			_singleFlipStatistic[i]++;
			return TakeResult();
		}

		for (size_t j = i + 1; j < k; j++)
		{
			int bitPosition2 = _unfrozenBits[j];
			if (_x[bitPosition2] == originalCodeword[bitPosition2])
				continue;

			_x[bitPosition2] =! _x[bitPosition2];
			DecodeFrom(bitPosition2);

			if (_x == _codeword) {
				_doubleFlipStatistic[i][j]++;
				return TakeResult();
			}

			_x[bitPosition2] = !_x[bitPosition2];
			DecodeFrom(bitPosition2);
		}

		_x[bitPosition] = !_x[bitPosition];
		DecodeFrom(bitPosition);
	}

	return TakeResult();
}