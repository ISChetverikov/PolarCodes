#include <functional>

#include "../include/Exceptions.h"
#include "../include/BaseDecoder.h"
#include "../include/Domain.h"

BaseDecoder::BaseDecoder(PolarCode * codePtr) {
	_codePtr = codePtr;
	_sigma = 0;
	

	_operationsCount = 0;
	_normalizerOperationCount = 0;
}

double BaseDecoder::GetOperationsCount() {
	if (_normalizerOperationCount == 0)
		return 0;

	return (double)_operationsCount / _normalizerOperationCount;
}

void BaseDecoder::ClearOperationsCount() {
	_operationsCount = 0;
	_normalizerOperationCount = 0;
}

void BaseDecoder::SetSigma(double sigma) {
	_sigma = sigma;
}

void BaseDecoder::SetDecoderAnswer(std::vector<int> codeword) {
	_codeword = codeword;
}

domain BaseDecoder::GetDomain() {
#ifdef DOMAIN_LLR
	return LLR;
#elif DOMAIN_P1
	return P1;
#endif
}

#ifdef DOMAIN_LLR

#ifdef MINSUM
double BaseDecoder::f(double x, double y) {
	double sign = 1.0;

	if (x < 0) {
		sign *= -1;
		x *= -1;

		// operations count
		_operationsCount += 2;
		///////////////////
	}
	if (y < 0) {
		sign *= -1;
		y *= -1;

		// operations count
		_operationsCount += 2;
		///////////////////
	}

	// operations count
	_operationsCount += 2;
	///////////////////

	return ((x < y) ? x : y) * sign;
}
#else
double BaseDecoder::f(double llr1, double llr2) {
	double prod = tanh(llr1 / 2) * tanh(llr2 / 2);
	double limit = 0.9999999999999999;

	if (prod > limit)
		prod = limit;
	if (prod < -limit)
		prod = -limit;

	// operations count
	_operationsCount += 9;
	///////////////////

	return 2 * atanh(prod);
}

#endif // MINSUM

double BaseDecoder::g(double x, double y, int b) {

	// operations count
	_operationsCount += 4;
	///////////////////

	return y + (1 - 2 * b) * x;
}

int BaseDecoder::L(double llr) {
	// operations count
	_operationsCount += 1;
	///////////////////

	return llr < 0;
}
#elif DOMAIN_P1

double BaseDecoder::f(double p1, double p2) {
	// operations count
	_operationsCount += 5;
	///////////////////

	return p1 * (1 - p2) + p2 * (1 - p1);
}

double BaseDecoder::g(double p1, double p2, int b) {
	double p1_b = f(b, p1);

	// operations count
	_operationsCount += 7;
	///////////////////

	if (p1_b == 0 && p2 == 1 || p2 == 0 && p1_b == 1)
		return 1 / 2;

	return p1_b * p2 / (p1_b * p2 + (1 - p1_b) * (1 - p2));
}

int BaseDecoder::L(double p1) {

	// operations count
	_operationsCount += 1;
	///////////////////

	return p1 >= (0.5);
}
#endif


std::string BaseDecoder::GetStatistic() {
	return "";
}
void BaseDecoder::ClearStatistic() {
	
}

std::string BaseDecoder::GetPathInfo() {
	return _pathTrace;
}
