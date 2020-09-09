#include <functional>

#include "../include/Exceptions.h"
#include "../include/BaseDecoder.h"
#include "../include/Domain.h"

BaseDecoder::BaseDecoder(PolarCode * codePtr, domain domain, bool isMinSum) {
	_codePtr = codePtr;
	_sigma = 0;
	_domain = domain;

	if (domain == LLR) {
		_g = [this](double llr1, double llr2, int b) {
			return this->g_Llr(llr1, llr2, b);
		};
		_L = [this](double llr) {
			return this->L_Llr(llr);
		};
		if (isMinSum)
			_f = [this](double llr1, double llr2) {
				return this->f_LlrMinSum(llr1, llr2);
			};
		else
			_f = [this](double llr1, double llr2) {
				return this->f_LlrTanh(llr1, llr2);
			};
	}
	else {
		_g = [this](double llr1, double llr2, int b) {
			return this->g_P1(llr1, llr2, b);
		};
		_L = [this](double llr) {
			return this->L_P1(llr);
		};
		_f = [this](double llr1, double llr2) {
			return this->f_P1(llr1, llr2);
		};
	}
	
}

double BaseDecoder::GetSigma() {
	return _sigma;
}

void BaseDecoder::SetSigma(double sigma) {
	_sigma = sigma;
}

domain BaseDecoder::GetDomain() {
	return _domain;
}

// LLR domain functions
double BaseDecoder::f_LlrMinSum(double x, double y) {
	double sign = 1.0;

	if (x < 0) {
		sign *= -1;
		x *= -1;
	}
	if (y < 0) {
		sign *= -1;
		y *= -1;
	}

	return ((x < y) ? x : y) * sign;
}

double BaseDecoder::f_LlrTanh(double llr1, double llr2) {
	double prod = tanh(llr1 / 2) * tanh(llr2 / 2);
	double limit = 0.9999999999999999;

	if (prod > limit)
		prod = limit;
	if (prod < -limit)
		prod = -limit;

	return 2 * atanh(prod);
}

double BaseDecoder::g_Llr(double x, double y, int b) {
	return y + (1 - 2 * b) * x;
}

int BaseDecoder::L_Llr(double llr) {
	return llr < 0;
}

// P1 domain functions
double BaseDecoder::f_P1(double p1, double p2) {
	return p1 * (1 - p2) + p2 * (1 - p1);
}

double BaseDecoder::g_P1(double p1, double p2, int b) {
	return g_est(f_P1(b, p1), p2);
}
double BaseDecoder::g_est(double p1, double p2) {
	if (p1 == 0 && p2 == 1 || p2 == 0 && p1 == 1)
		return 1 / 2;

	return p1 * p2 / (p1 * p2 + (1 - p1) * (1 - p2));
}

int BaseDecoder::L_P1(double p1) {
	return p1 >= (0.5);
}
