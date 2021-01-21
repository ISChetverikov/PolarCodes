#include "../include/ScListFanoDecoder.h"
#include "../include/GaussianApproximation.h"

ScListFanoDecoder::ScListFanoDecoder(PolarCode * code, double T, double delta, double approximationSnr, int L) : ScCrcAidedDecoder(code) {
	size_t n = _codePtr->N();

	_T = T;
	_delta = delta;
	_L = L;
	double approximationSigma = sqrt(pow(10, -1.0 / 10.0));
	GaussianApproximation ga(approximationSigma);
	_p = std::vector<double>(n, 0);
	for (size_t i = 0; i < n; i++)
	{
		_p[i] = ga.GetChannelErrorProbability(i + 1, n);
	}
}

void ScListFanoDecoder::UpdateT(double & T, double & delta, double & tau) {
	while (T + delta < tau)
		T += delta;
}

void ScListFanoDecoder::BackwardMove(std::vector<double> & beta, std::vector<bool> & gamma, int & i, double & T, double & delta, bool & B, int & j) {

	while (true) {
		double mu = 0;

		if (j <= -1)
			mu = -1000;

		if (j >= 1)
			mu = beta[j - 1];

		if (mu >= T) {
			j--;
			if (!gamma[j + 1])
			{
				B = true;
				return;
			}
		}
		else {
			T -= delta;
			B = false;
			return;
		}

	}
}

std::vector<int> ScListFanoDecoder::Decode(std::vector<double> beliefs) {
	size_t n = beliefs.size();
	size_t m = _codePtr->m();
	size_t k = _codePtr->k();
	for (size_t i = 0; i < n; i++)
	{
		_beliefTree[0][i] = beliefs[i];
	}

	int i = 0;
	int j = -1;
	bool B = false;
	std::vector<double> beta(k, 0.0); // only for frozen bits
	std::vector<double> metrics(n, 0); // for all bits
	std::vector<bool> gamma(k, 0);
	std::vector<int> A = _codePtr->UnfrozenBits(); // info set
	double T = _T;
	double delta = _delta;

	size_t firstInfoBit = 0;
	for (size_t i = 0; i < n; i++)
	{
		if (_maskWithCrc[i]) {
			firstInfoBit = i;
			break;
		}
	}
	while (i < n)
	{
		PassDown(i); // get p1 metric in _beliefTree[m][i]
		double p0 = 1 - _beliefTree[m][i];
		double p1 = _beliefTree[m][i];

		if (_maskWithCrc[i]) {
			double previous = (i == 0) ? 0 : metrics[i - 1];
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

					metrics[i] = max;
					beta[j + 1] = max;
					gamma[j + 1] = false;

					double mu = 0;
					if (j != -1)
						mu = beta[j];

					if (mu < T + delta)
						UpdateT(T, delta, beta[j + 1]);
					i++;
					j++;
				}
				else {
					if (min > T) {
						_x[i] = argmin;
						_uhatTree[m][i] = _x[i];
						PassUp(i);

						metrics[i] = min;
						beta[j + 1] = min;
						gamma[j + 1] = true;
						B = false;

						i++;
						j++;
					}
					else {
						if (j == -1) {
							T = T - delta;
							B = false;
						}
						else {
							BackwardMove(beta, gamma, i, T, delta, B, j);
							i = A[j + 1];
						}
					}
				}
			}
			else {
				if (j == -1)
					T = T - delta;

				else {
					BackwardMove(beta, gamma, i, T, delta, B, j);
					i = A[j + 1];
				}

			}
		}
		else {
			_x[i] = 0;
			_uhatTree[m][i] = _x[i];
			PassUp(i);

			double currentMetric = log(p0) - log(1 - _p[i]);
			// cumulative
			metrics[i] = currentMetric;
			if (i != 0)
				metrics[i] += metrics[i - 1];

			i++;
		}
	}

	return TakeResult();
}