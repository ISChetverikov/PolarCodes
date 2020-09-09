#pragma once

#include "ScCrcAidedDecoder.h"

class ScListDecoder : public ScCrcAidedDecoder {
private:
	int _L;

public:
	ScListDecoder(PolarCode * code, int L);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListDecoder() {};
};