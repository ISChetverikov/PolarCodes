#include "../include/PolarCode.h"
#include "../include/Exceptions.h"

PolarCode::PolarCode() {

}
PolarCode::PolarCode(int m, int k, std::vector<int> reliabilitySequence) {
	_m = m;
	_k = k;
	_N = 2 << (m - 1);
	size_t sequenceLength = reliabilitySequence.size();
	if (_N != sequenceLength)
		throw IncorrectSequenceSizeException("Sequence length does not match with code length");

	_bitsMask = std::vector<int>(_N, 0);
	for (size_t i = sequenceLength - k; i < sequenceLength; i++)
	{
		_bitsMask[reliabilitySequence[i]] = 1;
	}
}

size_t PolarCode::m() {
	return _m;
}
size_t PolarCode::N() {
	return _N;
}
size_t PolarCode::k() {
	return _k;
}
std::vector<int> PolarCode::BitsMask() {
	return _bitsMask;
}