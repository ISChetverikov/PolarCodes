#pragma once

#include "ScCrcAidedDecoder.h"

class ScFlipStatDecoder : public ScCrcAidedDecoder {
private:
	int _T;

public:
	ScFlipDecoder(PolarCode * code, int T);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScFlipDecoder() {};
};