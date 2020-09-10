#pragma once

#include "ScCrcAidedDecoder.h"

class ScListDecoder : public ScCrcAidedDecoder {
protected:
	int _L;
	std::vector<std::vector<std::vector<double>>> _beliefTrees;
	std::vector<std::vector<std::vector<int>>> _uhatTrees;
	std::vector<std::vector<int>> _candidates;
	std::vector<double> _metrics;

	void PassDownList(size_t iter);
	void PassUpList(size_t iter);
	void DecodeListInternal(std::vector<double> inLlr);
	double StepMetric(double belief, int decision);
public:
	ScListDecoder(PolarCode * code, int L);
	std::vector<int> Decode(std::vector<double> llr) override;
	~ScListDecoder() {};
};