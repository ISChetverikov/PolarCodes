
#include "../include/PolarCode.h"
#include "../include/ScDecoder.h"
#include "../include/ScListFanoDecoder.h"
#include "../include/Exceptions.h"
#include "../include/Domain.h"
#include "../include/GaussianApproximation.h"

#define DBL_MAX 1.7976931348623158e+308 
#define FROZEN_VALUE 0

ScListFanoDecoder::ScListFanoDecoder(PolarCode * code, double T, double delta, double approximationSnr, double L) : ScCrcAidedDecoder(code) {
	size_t n = _codePtr->N();
	size_t k = _codePtr->kExt();

	_T = T;
	_delta = delta;
	_L = L;

	double approximationSigma = sqrt(pow(10, -approximationSnr / 10.0) / 2);
	GaussianApproximation ga(approximationSigma);
	_p = std::vector<double>(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		_p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}

	// Inners states
	_beta = std::vector<double>(k, 0.0); // only for unfrozen bits
	_metrics = std::vector<double>(n, 0); // for all bits
	_gamma = std::vector<bool>(k, 0);
	_A = _codePtr->UnfrozenBitsWithCrc(); // info set
}

void ScListFanoDecoder::UpdateT(double & T, double & tau) {
	while (T + _delta < tau)
		T += _delta;
}

void ScListFanoDecoder::BackwardMove(double & T, bool & B, int & j, int rootIndex) {

	while (true) {
		double mu = 0;

		if (j <= rootIndex - 1)
			mu = -1000;

		if (j >= 1)
			mu = _beta[j - 1];

		if (mu >= T) {
			j--;
			if (!_gamma[j + 1])
			{
				B = true;
				return;
			}
		}
		else {
			T -= _delta;
			B = false;
			return;
		}

	}
}

void ScListFanoDecoder::DecodeFrom(int rootIndex) {
	size_t n = _codePtr->N();
	size_t m = _codePtr->m();
	size_t k = _codePtr->kExt();
	
	int j = rootIndex;
	int i = 0;
	if (j >= 0)
		i = _A[j + 1];

	bool B = false;
	double T = _T;

	while (i < n)
	{
		PassDown(i); // get p1 metric in _beliefTree[m][i]
		double p0 = 1 - _beliefTree[m][i];
		double p1 = _beliefTree[m][i];

		if (_maskWithCrc[i]) {
			double previous = (i == 0) ? 0 : _metrics[i - 1];
			double m0 = previous + log(p0 / (1 - _p[i]));
			double m1 = previous + log(p1 / (1 - _p[i]));

			double max = (m1 > m0) ? m1 : m0;
			int argmax = (m1 > m0) ? 1 : 0;

			double min = (m1 > m0) ? m0 : m1;
			int argmin = (m1 > m0) ? 0 : 1;

			if (max > T) {
				if (!B) {
					_x[i] = argmax;
					_uhatTree[m][i] = _x[i];
					PassUp(i);

					_metrics[i] = max;
					_beta[j + 1] = max;
					_gamma[j + 1] = false;

					double mu = 0;
					if (j != -1)
						mu = _beta[j];

					if (mu < T + _delta)
						UpdateT(T, _beta[j + 1]);
					i++;
					j++;
				}
				else {
					if (min > T) {
						_x[i] = argmin;
						_uhatTree[m][i] = _x[i];
						PassUp(i);

						_metrics[i] = min;
						_beta[j + 1] = min;
						_gamma[j + 1] = true;
						B = false;

						i++;
						j++;
					}
					else {
						if (j == -1) {
							T = T - _delta;
							B = false;
						}
						else {
							BackwardMove(T, B, j, rootIndex);
							i = _A[j + 1];
						}
					}
				}
			}
			else {
				if (j == -1)
					T = T - _delta;

				else {
					BackwardMove(T, B, j, rootIndex);
					i = _A[j + 1];
				}

			}
		}
		else {
			_x[i] = FROZEN_VALUE;
			_uhatTree[m][i] = _x[i];
			PassUp(i);

			double currentMetric = log(p0) - log(1 - _p[i]);
			// cumulative
			_metrics[i] = currentMetric;
			if (i != 0)
				_metrics[i] += _metrics[i - 1];

			i++;
		}
	}
}

std::vector<int> ScListFanoDecoder::Decode(std::vector<double> beliefs) {
	size_t n = beliefs.size();
	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = beliefs[i];
	}

	int rootIndex = -1; // j started position
	DecodeFrom(rootIndex);

	return TakeResult();
}