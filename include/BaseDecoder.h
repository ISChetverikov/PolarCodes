#pragma once

#include <vector>
#include <functional>

#include "PolarCode.h"
#include "Domain.h"

class BaseDecoder {

protected:
	PolarCode * _codePtr;
	double _sigma;
	
	double f(double p1, double p2);
	double g(double p1, double p2, int b);
	int L(double p1);

public:
	BaseDecoder(PolarCode * codePtr);
	
	virtual domain GetDomain();
	virtual void SetSigma(double sigma);
	double GetSigma();
	virtual std::vector<int> Decode(std::vector<double> llr) = 0;
	~BaseDecoder() {};
};