#include <mpi.h>
#include <iostream>
#include <thread>
#include "StatisticalSequence.h"

using std::cout;
using std::string;
int main(int argc, char * argv[]) {

	int ProcNum, ProcRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	int maxTestsCount = 10000;
	int maxRejectionsCount = 50;
	
	//auto folder = "C:\\Users\\ische\\source\\repos\\PolarCodes\\polar_sequences\\Stat\\OpenMPI\\";
	string folder = "";
	double snr = 0.0;
	int m = 0;
	int k = 0;
	// double snr = PickUpSnr(m, k);

	int isChecked = 1;
	if (ProcRank == 0) {
		if (argc != 4) {
			cout << "Usage: " << argv[0] << " <ouptut folder> <m> <k>" << std::endl;
			isChecked = 0;
		}
		else {
			m = std::stoi(argv[2]);
			k = std::stoi(argv[3]);

			folder = argv[1];
		}
	}
	MPI_Bcast(&isChecked, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (isChecked == 0) {
		MPI_Finalize();
		return 1;
	}
	MPI_Bcast(&snr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);


	BuiltSequenceStatistically(folder, m, k, maxTestsCount, maxRejectionsCount, snr, ProcRank, ProcNum);
	MPI_Finalize();

	return 0;

}