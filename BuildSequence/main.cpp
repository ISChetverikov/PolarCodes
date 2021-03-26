#include <mpi.h>
#include <iostream>
#include "StatisticalSequence.h"

using std::cout;

int main(int argc, char * argv[]) {

	int ProcNum, ProcRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	auto folder = "C:\\Users\\ische\\source\\repos\\PolarCodes\\polar_sequences\\Stat\\";
	int maxTestsCount = 1000;
	int maxRejectionsCount = 10;
	int m = 8;
	int k = 128;
	double snr = 0.5;
	//double snr = PickUpSnr(m, k);



	if (ProcRank == 0) {
		if (argc != 3) {
			cout << "Usage: " << argv[0] << " <input file> <output file>" << std::endl;
			isChecked = 0;
		}
	}

	BuiltSequenceStatistically(folder, m, k, maxTestsCount, maxRejectionsCount, snr);
	
	MPI_Finalize();
	return 0;

}