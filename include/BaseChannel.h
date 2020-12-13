#pragma once

#include <vector>

class BaseChannel {
protected:
	double _sigma = 1.0;

public:
	BaseChannel();
	virtual ~BaseChannel() {};
	virtual void SetSigma(double sigma);
	virtual std::vector<double> Pass(std::vector<int> codeword) = 0;
};