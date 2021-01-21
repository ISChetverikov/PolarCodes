#include "../include/ScListFanoDecoder.h"
#include "../include/GaussianApproximation.h"

ScListFanoDecoder::ScListFanoDecoder(PolarCode * code, double T, double delta, double approximationSnr, int L) : ScCrcAidedDecoder(code) {
	size_t n = _codePtr->N();
	size_t k = _codePtr->kExt();

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

	_betaJ = std::vector<double>(k, 0.0); // only for frozen bits
	_betaI = std::vector<double>(n, 0); // for all bits
	_gamma = std::vector<bool>(k, 0);
}

FanoState ScListFanoDecoder::SaveState() {
	FanoState state = {};

	return state;
}

void ScListFanoDecoder::LoadState(FanoState state) {

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

	_i = 0;
	_j = -1;
	_B = false;
	
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
	while (_i < n)
	{
		PassDown(_i); // get p1 metric in _beliefTree[m][i]
		double p0 = 1 - _beliefTree[m][_i];
		double p1 = _beliefTree[m][_i];

		if (_maskWithCrc[_i]) {
			double previous = (_i == 0) ? 0 : _betaI[_i - 1];
			double m0 = previous + log(p0 / (1 - _p[_i]));
			double m1 = previous + log(p1 / (1 - _p[_i]));

			double max = (m1 > m0) ? m1 : m0;
			int argmax = (m1 > m0) ? 1 : 0;

			double min = (m1 > m0) ? m0 : m1;
			int argmin = (m1 > m0) ? 0 : 1;

			if (max > T) {
				if (!_B) {
					_x[_i] = argmax;
					_uhatTree[m][_i] = _x[_i];
					PassUp(_i);

					_betaI[_i] = max;
					_betaJ[_j + 1] = max;
					_gamma[_j + 1] = false;

					double mu = 0;
					if (_j != -1)
						mu = _betaJ[_j];

					if (mu < T + delta)
						UpdateT(T, delta, _betaJ[_j + 1]);
					_i++;
					_j++;
				}
				else {
					if (min > T) {
						_x[_i] = argmin;
						_uhatTree[m][_i] = _x[_i];
						PassUp(_i);

						_betaI[_i] = min;
						_betaJ[_j + 1] = min;
						_gamma[_j + 1] = true;
						_B = false;

						_i++;
						_j++;
					}
					else {
						if (_j == -1) {
							T = T - delta;
							_B = false;
						}
						else {
							BackwardMove(_betaJ, _gamma, _i, T, delta, _B, _j);
							_i = A[_j + 1];
						}
					}
				}
			}
			else {
				if (_j == -1)
					T = T - delta;

				else {
					BackwardMove(_betaJ, _gamma, _i, T, delta, _B, _j);
					_i = A[_j + 1];
				}

			}
		}
		else {
			_x[_i] = 0;
			_uhatTree[m][_i] = _x[_i];
			PassUp(_i);

			double currentMetric = log(p0) - log(1 - _p[_i]);
			// cumulative
			_betaI[_i] = currentMetric;
			if (_i != 0)
				_betaI[_i] += _betaI[_i - 1];

			_i++;
		}
	}

	return TakeResult();
}