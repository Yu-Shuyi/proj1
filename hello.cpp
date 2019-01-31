#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int dboard(int N) {
	int M = 0;
	double r = 1;
	for (int i = 0; i < N; i++){
		double a = r * r * (rand() % 1000 / 999.0);
		double theta = rand();
		double x = sqrt(a) * cos(theta);
		double y = sqrt(a) * sin(theta);
		if (x <= r / sqrt(2) && y <= r / sqrt(2) && x >= -r / sqrt(2) && y >= -r / sqrt(2)){
			M++;
		}
	}
	return M;
}

int main(int argc, char *argv[]) {
	// set up MPI
	MPI_Init(&argc, &argv);
	double t_start = MPI_Wtime();
	
	// get communicator size and my rank
	MPI_Comm comm = MPI_COMM_WORLD;
	int P, rank;
	MPI_Comm_size(comm, &P);
	MPI_Comm_rank(comm, &rank);
	
	// !!!did not check input validity!!!
	int N, R;
	if (rank == 0) sscanf(argv[1], "%d", &N);
	sscanf(argv[2], "%d", &R);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	srand(rank);
	int m = 0, n = N / P;
	if (rank < N % P) n += 1;
	for (int i = 0; i < R; i++){
		m += dboard(n);
	}
	int sum_m = 0;
	MPI_Reduce(&m, &sum_m, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	double pi;
	if (rank == 0) {
		pi = 2.0 * N *R / sum_m;
	}
	
	// finalize MPI
	MPI_Barrier(MPI_COMM_WORLD);
	double t_end = MPI_Wtime();
	if (rank == 0) printf("N = %d, R = %d, P = %d, PI = %f\nTime = %f\n", N, R, P, pi, t_end - t_start);
	return MPI_Finalize();
}
