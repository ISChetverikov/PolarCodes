#pragma once

#include "ScCrcAidedDecoder.h"

struct FanoState {
	int i;
	int j;

	std::vector<std::vector<double>> beliefTree;
	std::vector<std::vector<int>> uhatTree;
	
	std::vector<double> beta;
	std::vector<double> metrics;
	std::vector<bool> gamma;
};

class ScListFanoDecoder : public ScCrcAidedDecoder {

private:
	double _T;
	double _delta;
	int _L;

	std::vector<double> _p; // channel error probabilities 

	// list stuff
	std::vector<FanoState> states;
	std::vector<double> metrics;

	void UpdateT(double & T, double & delta, double & tau);
	void BackwardMove(std::vector<double> & beta, std::vector<bool> & gamma, int & i, double & T, double & delta, bool & B, int & j);

public:
	ScListFanoDecoder(PolarCode * code, double T, double delta, double approximationSnr, int L);
	// All fano decoders works only with P1 domain
	std::vector<int> Decode(std::vector<double> p1) override;
	~ScListFanoDecoder() {};
};