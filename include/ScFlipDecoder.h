#pragma once

#include "ScCrcAidedDecoder.h"

class ScFlipDecoder : public ScCrcAidedDecoder {
private:
	int _T;

protected:
	
	std::vector<int> GetSmallestLlrsIndices(std::vector<double> llrs, int count);

public:
	ScFlipDecoder(PolarCode * code, domain domain, bool isMinSum, int T);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScFlipDecoder() {};
};