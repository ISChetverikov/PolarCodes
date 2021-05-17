#pragma once

#include "BaseDecoder.h"

using std::vector;

class ScStackOptimizedDecoder : public BaseDecoder {

protected:

	size_t _L = 0;
	size_t _D = 0;

	// _n = 2 ^ _m
	size_t _m = 0;
	size_t _n = 0;

	vector<vector<double>> _alpha;
	vector<vector<vector<int>>> _beta;

	vector<int> _mask;

	double f(double left, double right);
	double g(double left, double right, int left_hard);
	int HD(double llr);

	void recursively_calc_alpha(size_t lambda, size_t phi);
	void recursively_calc_beta(size_t lambda, size_t phi);

public:
	ScStackOptimizedDecoder(PolarCode * codePtr, int L, int D);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScStackOptimizedDecoder() {};
};