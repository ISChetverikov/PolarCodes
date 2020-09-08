#pragma once

#include <vector>
#include <functional>

#include "PolarCode.h"
#include "Domain.h"

class BaseDecoder {

protected:
	PolarCode * _codePtr;
	double _sigma;
	domain _domain;
	std::function<double(double,double)> _f;
	std::function<double(double,double,int)> _g;
	std::function<int(double)> _L;

	// LLR funcitons (f - left estimation, g - right estimation, L - soft to hard desicion)
	double f_LlrMinSum(double x, double y);
	double f_LlrTanh(double llr1, double llr2);
	double g_Llr(double x, double y, int b);
	int L_Llr(double llr);

	// P1 functions ( + g_est - bit-node operation separately)
	double f_P1(double p1, double p2);
	double g_P1(double p1, double p2, int b);
	double g_est(double p1, double p2);
	int L_P1(double p1);

public:
	BaseDecoder(PolarCode * codePtr, domain domain, bool isMinSum);
	
	virtual domain GetDomain();
	virtual void SetSigma(double sigma);
	double GetSigma();
	virtual std::vector<int> Decode(std::vector<double> llr) = 0;
	~BaseDecoder() {};
};