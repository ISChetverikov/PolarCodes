#pragma once

#include <deque>

#include "BaseDecoder.h"
#include "ScOptimized.h"

using std::vector;
using std::deque;
using std::pair;

class ScCreeperDecoder : public ScOptimizedDecoder {

protected:;
	double _delta;

	// nodes as pair of (bit number, bit value)
	deque<pair<int, int>> _NP;
	deque<double> _TP;
	deque<bool> _F;

	vector<double> _metric;

	void recursively_calc_alpha_creeper(size_t lambda, size_t phi, bool isPathSwitched);
	double calculate_step_metric(double newLlr, int decision);
	double Q(double x, double T);
public:
	ScCreeperDecoder(PolarCode * codePtr, double delta);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScCreeperDecoder() {};
};