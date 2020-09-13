#include <sstream>

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

	ss << "Single Flip:\n";
	for (size_t i = 0; i < _codePtr->k(); i++)
	{
		ss << "(" << _unfrozenBits[i] << "): " << _singleFlipStatistic[i] << "\n";
	}

	for (size_t i = 0; i < _codePtr->k(); i++)
	{
		for (size_t j = 0; j < _codePtr->k(); j++)
		{
			ss << "(" << _unfrozenBits[i] << ", " << _unfrozenBits[j] << "): " << _doubleFlipStatistic[i][j] << "\n";
		}
	}

}

std::vector<int>  ScFlipStatDecoder::Decode(std::vector<double> beliefs) {
	DecodeInternal(beliefs);
	
	auto originalCodeword = _codeword;

	if (_x == originalCodeword)
		return TakeResult();

	size_t k = _codePtr->k();
	for (size_t i = 0; i < k; i++)
	{

	}
}