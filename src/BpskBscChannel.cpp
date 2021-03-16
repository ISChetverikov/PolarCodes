#include <iostream>
#include <random>
#include "../include/BpskBscChannel.h"
#include "../include/CommonTransformations.h"

BpskBscChannel::BpskBscChannel(double coderate) : BaseChannel() {
	std::bernoulli_distribution bernoulli_dist(0.5);
	_bernoulli_dist = bernoulli_dist;
	
	_coderate = coderate;
	_fixedLlr = 1.0;
}


void BpskBscChannel::SetSnr(double snr) {
	BaseChannel::SetSnr(snr);

	//double p = ebnoToPErr(snrToEbN0(_snr, _coderate));
	double p = ebnoToPErr(snrToSigma(snr));

	_fixedLlr = log((1 - p) / p);
	std::bernoulli_distribution b(p);
	_bernoulli_dist = b;
}


std::vector<double> BpskBscChannel::Pass(std::vector<int> codeword) {
	size_t n = codeword.size();
	std::vector<double> output(n, 0);

	for (size_t i = 0; i < n; i++) {

		if (_bernoulli_dist(_randomDevice)) {
			output[i] = !codeword[i];
		}
		else
		{
			output[i] = codeword[i];
		}

		if (output[i] == 1) {
			output[i] = -_fixedLlr;
		}
		else {
			output[i] = _fixedLlr;
		}
	}

#ifdef DOMAIN_P1

	for (size_t i = 0; i < n; i++)
	{
		output[i] = LlrToP12(output[i]);
	}

#endif // DOMAIN_P1

	return output;

}