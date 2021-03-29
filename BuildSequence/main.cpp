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

	cout << "Start sleep...\n" << std::flush;
	//std::this_thread::sleep_for(std::chrono::milliseconds(20000));
	cout << "End sleep...\n";

	int maxTestsCount = 1000;
	int maxRejectionsCount = 10;
	
	//auto folder = "C:\\Users\\ische\\source\\repos\\PolarCodes\\polar_sequences\\Stat\\OpenMPI\\";
	string folder = "";
	double snr = 0.0;
	int m = 5;
	int k = 16;
	//double snr = PickUpSnr(m, k);

	int isChecked = 1;
	if (ProcRank == 0) {
		if (argc != 5) {
			cout << "Usage: " << argv[0] << " <folder> <snr> <m> <k>" << std::endl;
			isChecked = 0;
		}
		else {
			snr = std::stod(argv[2]);
			m = std::stoi(argv[3]);
			k = std::stoi(argv[4]);

			folder = argv[1];
		}
	}
	MPI_Bcast(&isChecked, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "Rank: " << ProcRank << " isChecked: " << isChecked << std::endl;
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