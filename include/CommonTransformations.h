#pragma once


double snrToSigma(double snr);

double snrToEbN0(double snr, double coderate);

double ebnoToPErr(double sigma);

double LlrToP1(double llr);
