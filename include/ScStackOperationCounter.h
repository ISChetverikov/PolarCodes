#pragma once
#include <stack>

#include "CRC.h"
#include "BaseDecoder.h"
#include "ScOptimized.h"

using std::vector;
using std::stack;

class ScStackOperationCounter : public ScOptimizedDecoder {

protected:

	size_t _L = 0;
	size_t _D = 0;

	size_t _kWithCrc;
	CRC * _crcPtr;

	vector<double> _path_metrics;
	vector<size_t> _path_lengths;
	vector<size_t> _path_lengths_hits;
	vector<vector<int>> _paths;
	stack<size_t> _inactive_path_indices;
	vector<bool> _is_paths_active;

	double f(double left, double right);
	double g(double left, double right, int left_hard);
	int HD(double llr);

	void recursively_calc_alpha_stack(size_t lambda, size_t phi, bool isPathSwitched);
	void recursively_calc_beta(size_t lambda, size_t phi);

	double calculate_step_metric(double newLlr, int decision);
	bool IsCrcPassed(vector<int> & word);

	void KillPath(size_t index, size_t & T);
	size_t ClonePath(size_t index, size_t & T);
	void ExtendPath(size_t index, int bit, double pathMetric);
	void LoadPath(size_t index);

public:
	ScStackOperationCounter(PolarCode * codePtr, int L, int D);

	std::vector<int> Decode(std::vector<double> llr) override;

	~ScStackOperationCounter() {};
};