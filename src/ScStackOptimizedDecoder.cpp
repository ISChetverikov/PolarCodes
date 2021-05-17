#include "../include/PolarCode.h"
#include "../include/ScStackOptimizedDecoder.h"

#define FROZEN_VALUE 0
#define MAX_METRIC 1000000.0

ScStackOptimizedDecoder::ScStackOptimizedDecoder(PolarCode * codePtr, int L, int D) : ScOptimizedDecoder(codePtr) {
	_L = L;
	_D = D;

	_kWithCrc = _codePtr->kExt();
	_crcPtr = new CRC(_codePtr->CrcPoly());

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

void ScStackOptimizedDecoder::recursively_calc_alpha_stack(size_t lambda, size_t phi, bool isPathSwitched) {
	_operationsCount += 1;

	if (lambda == _m)
		return;

	_operationsCount += 4;

	size_t lambda_big = 1 << lambda;
	size_t lambda_next = lambda + 1;
	size_t isPhiOdd = phi % 2;

	if (!isPhiOdd || isPathSwitched)
		recursively_calc_alpha(lambda_next, phi >> 1);

	if (isPhiOdd) {
		for (size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = g(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big], _beta[!isPhiOdd][lambda][i]);
			_operationsCount += 3;
		}
	}
	else {

		_operationsCount += 1;

		for (size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = f(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big]);

			_operationsCount += 3;
		}
	}

	return;
}

double ScStackOptimizedDecoder::calculate_step_metric(double newLlr, int decision) {
	// economize operations
	return (newLlr < 0 && decision == 0 || newLlr > 0 && decision == 1) ? fabs(newLlr) : 0;
}

// Here word without frozen bits, all crc bits are located in the begging of the word
bool ScStackOptimizedDecoder::IsCrcPassed(vector<int> & codeword) {
	size_t deg = _codePtr->CrcDeg();

	vector<int> word(_k, 0);
	vector<int> crc(deg, 0);

	for (size_t i = 0; i < _k; i++)
	{
		word[i] = codeword[deg + i];
	}
	for (size_t i = 0; i < deg; i++)
	{
		crc[i] = codeword[i];
	}

	_crcPtr->ClearOperationsCount();
	auto crcReal = _crcPtr->Calculate(word);

	// operations count
	_operationsCount += 1 + _crcPtr->GetLastOperationsCount();
	///////////////////

	return crc == crcReal;
}


std::vector<int> ScStackOptimizedDecoder::Decode(std::vector<double> llr) {
	vector<int> result = vector<int>(_k, 0);

	// variables
	size_t T = 0;
	size_t min_index = 0;
	size_t max_index = 0;
	bool isPathSwitched = false;
	double pm0 = 0.0;
	double pm1 = 0.0;
	size_t phi = 0;
	bool isPhiEven = false;

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
		recursively_calc_alpha_stack(0, _path_lengths[min_index], isPathSwitched);

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

		// Update length info ?????
		if (poped elemtns count length <= L)
		// Check CRC and stack is empty
		if (path == n and IsCrcPassed() and stack is empty)
			break;

		if (isPathSwitched)
			LoadPath();

		phi = _path_lengths[min_index] - 1;
		_beta[phi % 2][0][0] = _paths[min_index][phi];
		if (phi % 2 == 1)
			recursively_calc_beta(0, phi);
	}

	std::vector<int> codewordBits = _codePtr->UnfrozenBits();
	for (size_t i = 0; i < codewordBits.size(); i++)
	{
		result[i] = _paths[min_index][codewordBits[i]];

		// operations count
		_operationsCount += 2;
		///////////////////
	}

	_normalizerOperationCount++;

	return result;
}