

#include "ScCrcAidedDecoder.h"
#pragma once

class SCIvanDecoder : public ScCrcAidedDecoder {
protected:
	int _L;
	int _k;
	int _find;
	std::vector<std::vector<std::vector<double>>> _beliefTrees;
	std::vector<std::vector<std::vector<int>>> _uhatTrees;
	std::vector<std::vector<int>> _candidates;
	std::vector<std::vector<int>> _meta_candidates;
	std::vector<double> _metrics;
	std::vector <double> _meta_metrics;

	// optimization allocation
	std::vector<bool> _areTakenZero;
	std::vector<bool> _areTakenOne;

	void PassDownList(size_t iter);
	void PassUpList(size_t iter);
	void DecodeListInternal(std::vector<double> inLlr);
	void FillListMask(size_t iter);
	double StepMetric(double belief, int decision);
	// void ChangeFind(int val);
	std::vector<int> TakeListResult();
	std::vector<int> TakeListResultFinal();

public:
	SCIvanDecoder(PolarCode * code, int L, int k);
	std::vector<int> Decode(std::vector<double> llr) override;
	~SCIvanDecoder() {};
};