

#include "StatisticalSequence.h"

int main() {

	auto folder = "C:\\Users\\ische\\source\\repos\\PolarCodes\\polar_sequences\\Stat\\";
	int maxTestsCount = 1000;
	int maxRejectionsCount = 10;
	int m = 8;
	int k = 128;
	double snr = 0.5;
	//double snr = PickUpSnr(m, k);
	BuiltSequenceStatistically(folder, m, k, maxTestsCount, maxRejectionsCount, snr);

	return 0;

}