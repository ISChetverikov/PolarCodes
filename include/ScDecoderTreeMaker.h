#pragma once

#include "ScDecoder.h"

class ScDecoderTreeMaker : public ScDecoder {

private:
	std::vector<double> _p;

public:
	ScDecoderTreeMaker(PolarCode * code, double approximationSigma);
	std::vector<int> Decode(std::vector<double> llr) override;

	~ScDecoderTreeMaker() {};
};