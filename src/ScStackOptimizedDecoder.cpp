#include "../include/PolarCode.h"
#include "../include/ScStackOptimizedDecoder.h"

#define FROZEN_VALUE 0
#define MAX_METRIC 1000000.0

ScStackOptimizedDecoder::ScStackOptimizedDecoder(PolarCode * codePtr, int L, int D) : ScOptimizedDecoder(codePtr) {
	_L = L;
	_D = D;

	_kWithCrc = _codePtr->kExt();

	_path_metrics = vector<double>(_D, 0.0);
	_path_lengths = vector<size_t>(_D, 0);
	_path_lengths_hits = vector<size_t>(_D, 0);
	//_inactive_path_indices = vector<size_t>(_D, 0);
	_is_paths_active = vector<bool>(_D, false);

	vector<int> temp = vector<int>(_kWithCrc, 0);
	for (size_t i = 0; i < _D; i++)
	{
		_paths.push_back(temp);
	}
}

double ScStackOptimizedDecoder::calculate_step_metric(double newLlr, int decision) {
	// economize operations
	return (newLlr < 0 && decision == 0 || newLlr > 0 && decision == 1) ? fabs(newLlr) : 0;
}

std::vector<int> ScStackOptimizedDecoder::Decode(std::vector<double> llr) {

	// variables
	size_t T = 0;
	size_t min_index = 0;
	size_t max_index = 0;
	bool isPathSwitched = false;
	double pm0 = 0.0;
	double pm1 = 0.0;

	// init structures
	for (size_t i = 0; i < _D; i++)
	{
		_inactive_path_indices.push(i);
		_is_paths_active[i] = false;
	}
	for (size_t i = 0; i < _n; i++)
	{
		_path_lengths_hits[i] = 0;
	}

	// set initial path
	min_index = _inactive_path_indices.top();
	_inactive_path_indices.pop();
	_is_paths_active[min_index] = true;
	_path_metrics[min_index] = 0.0;
	_path_lengths[min_index] = 0;
	T = 1;

	// load input llrs
	for (size_t phi = 0; phi < _n; phi++)
	{
		_alpha[_m][phi] = llr[phi];
	}

	while (true) {
		recursively_calc_alpha(0, _path_lengths[min_index]);

		pm0 = _path_metrics[min_index] + calculate_step_metric(_alpha[0][0], 0);
		pm1 = _path_metrics[min_index] + calculate_step_metric(_alpha[0][0], 1);

		if (_mask[_path_lengths[min_index]]) {
			if (T == _D && _path_metrics[max_index] > (pm0 > pm1 ? pm0 : pm1)) {
				KillPath(max_index);
			}
			
			if (pm0 < pm1) {
				if (T < _D) {
					max_index = ClonePath(min_index);
					ExtendPath(max_index, 1, pm1);
				}

				ExtendPath(min_index, 0, pm0);
			}
			else {
				if (T < _D) {
					max_index = ClonePath(min_index);
					ExtendPath(max_index, 0, pm0);
				}

				ExtendPath(min_index, 1, pm1);
			}

		}
		else {
			ExtendPath(min_index, FROZEN_VALUE, pm0);
		}

		// loop for new indices
		size_t new_min_index = 0;
		size_t new_max_index = 0;
		double min_metric = MAX_METRIC;
		double max_metric = -MAX_METRIC;
		
		for (size_t i = 0; i < _D; i++)
		{
			if (_path_metrics[i] < min_metric) {
				min_metric = _path_metrics[i];
				new_min_index = i;
			}
			if (_path_metrics[i] > max_metric) {
				max_metric = _path_metrics[i];
				new_max_index = i;
			}
		}

		isPathSwitched = min_index == new_min_index;
		min_index = new_min_index;
		max_index = new_max_index;
	}

	// Update length info
	// Check CRC and stack is empty

}