#pragma once
#include <vector>
class PolarCode {

protected:
	
	size_t _m = 0;
	size_t _N = 0;
	size_t _k = 0;
	std::vector<int> _bitsMask;
	std::vector<int> _unfrozenBits;

public:
	PolarCode();
	PolarCode(int m, int k, std::vector<int> reliabilitySequence);

	size_t m();
	size_t N();
	size_t k();
	std::vector<int> BitsMask();
	std::vector<int> UnfrozenBits();

	~PolarCode() {};
};