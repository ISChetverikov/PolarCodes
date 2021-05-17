
#include "../include/PolarCode.h"
#include "../include/ScOptimized.h"

#define FROZEN_VALUE 0

ScOptimizedDecoder::ScOptimizedDecoder(PolarCode * codePtr) : BaseDecoder(codePtr) {
	_m = _codePtr->m();
	_n = _codePtr->N();
	_k = _codePtr->k();

	vector<vector<int>> _beta_temp;

	for (size_t i = 0; i <= _m; i++)
	{
		_alpha.push_back(vector<double>(1 << i, 0.0));
		_beta_temp.push_back(vector<int>(1 << i, 0.0));
	}

	_beta.push_back(_beta_temp);
	_beta.push_back(_beta_temp);

	_mask = _codePtr->BitsMaskWithCrc();
}

double ScOptimizedDecoder::f(double left, double right) {
	double sign = 1;

	if (left < 0) {
		left = -left;
		sign = -sign;

		_operationsCount += 2;
	}
	if (right < 0) {
		right = -right;
		sign = -sign;

		_operationsCount += 2;
	}

	_operationsCount += 4;

	return ((left < right) ? left : right) * sign;
}

double ScOptimizedDecoder::g(double left, double right, int left_hard) {
	_operationsCount += 3;

	return right + (left_hard == 1 ? -1 : 1) * left;
}

int ScOptimizedDecoder::HD(double llr) {
	_operationsCount += 1;

	return llr < 0;
}

void ScOptimizedDecoder::recursively_calc_alpha(size_t lambda, size_t phi) {
	_operationsCount += 1;

	if (lambda == _m)
		return;
	
	_operationsCount += 4;

	size_t lambda_big = 1 << lambda;
	size_t lambda_next = lambda + 1;
	size_t isPhiOdd = phi % 2;

	if (isPhiOdd) {
		for (size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = g(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big], _beta[!isPhiOdd][lambda][i]);
			_operationsCount += 3;
		}
	}
	else {
		recursively_calc_alpha(lambda_next, phi >> 1);

		_operationsCount += 1;

		for (size_t i = 0; i < lambda_big; i++) {
			_alpha[lambda][i] = f(_alpha[lambda_next][i], _alpha[lambda_next][i + lambda_big]);

			_operationsCount += 3;
		}
	}
		
	return;
}

void ScOptimizedDecoder::recursively_calc_beta(size_t lambda, size_t phi) {

	_operationsCount += 2;

	size_t isPhiOdd = phi % 2;
	size_t isPhiEven = !isPhiOdd;
	if (isPhiEven)
		return;

	_operationsCount += 5;

	size_t lambda_big = 1 << lambda;
	size_t lambda_next = lambda + 1;
	size_t isHalfPhiOdd = (phi >> 1) % 2;

	for (size_t i = 0; i < lambda_big; i++)
	{
		_beta[isHalfPhiOdd][lambda_next][i] = _beta[isPhiEven][lambda][i] ^ _beta[isPhiOdd][lambda][i];
		_beta[isHalfPhiOdd][lambda_next][lambda_big + i] = _beta[isPhiOdd][lambda][i];

		_operationsCount += 4;
	}

	recursively_calc_beta(lambda_next, phi >> 1);
}

std::vector<int> ScOptimizedDecoder::Decode(std::vector<double> llr) {
	std::vector<int> result(_codePtr->k(), 0);
	size_t i = 0;

	for (size_t phi = 0; phi < _n; phi++)
	{
		_alpha[_m][phi] = llr[phi];

		_operationsCount += 2;
	}

	for (size_t phi = 0; phi < _n; phi++)
	{
		recursively_calc_alpha(0, phi);

		if (_mask[phi]) {
			int hard_desicion = HD(_alpha[0][0]);

			_beta[phi % 2][0][0] = hard_desicion;
			result[i++] = hard_desicion;

			_operationsCount += 2;
		}
		else {
			_beta[phi % 2][0][0] = FROZEN_VALUE;

			_operationsCount += 1;
		}

		recursively_calc_beta(0, phi);

		_operationsCount += 2;
	}

	_normalizerOperationCount++;

	return result;
}