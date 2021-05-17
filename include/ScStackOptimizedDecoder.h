#pragma once

#include <stack>

#include "BaseDecoder.h"
#include "ScOptimized.h"

using std::vector;
using std::stack;

class ScStackOptimizedDecoder : public ScOptimizedDecoder {

protected:

	size_t _L = 0;
	size_t _D = 0;

	size_t _kWithCrc;

	vector<double> _path_metrics;
	vector<size_t> _path_lengths;
	vector<size_t> _path_lengths_hits;
	vector<vector<int>> _paths;
	stack<size_t> _inactive_path_indices;
	vector<bool> _is_paths_active;

	double calculate_step_metric(double newLlr, int decision);

public:
	ScStackOptimizedDecoder(PolarCode * codePtr, int L, int D);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScStackOptimizedDecoder() {};
};