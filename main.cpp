#include <iostream>
#include <mpi.h>
#include <cmath>

#include <iostream>
#include <cstring>

#define MASTER 0

int N;
int P;
int Q;
int N_over_Q;

int* completeMatrix = nullptr;
int** subMatrix;

void min_plus_matrix_multiply(int* own, const int* other);

int main(int argc, char *argv[]) {
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MASTER) {

        std::cin >> N;

        Q = static_cast<int>(sqrt(P));
        if ((Q * Q) != P || N % Q != 0) {
            std::cout << "Incompatible values for P and N!" << std::endl;
            return EXIT_FAILURE;
        }

        N_over_Q = N / Q;
        completeMatrix = new int[N * N];
        std::memset(completeMatrix, 0, N*N);
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                int i, j, k;
                i = y/N_over_Q;
                j = x/N_over_Q;
                k = (y % N_over_Q)*N_over_Q + (x % N_over_Q);

                std::cin >> completeMatrix[i*N_over_Q*N + j*N_over_Q*N_over_Q + k];
            }
        }

        for (int i = 0; i < P; ++i) {
            std::cout << i << ":\t";
            for (int j = 0; j < N_over_Q*N_over_Q; ++j) {
                std::cout << completeMatrix[i*N_over_Q*N_over_Q + j] << " ";
            }
            std::cout << std::endl;
        }


    }
    // broadcast N/Q
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Q, 1, MPI_INT, 0, MPI_COMM_WORLD);
    N_over_Q = N / Q;

    // create cartesian grid communicator
    int grid_rank, reorder, dims[2], periods[2], coords[2];
    MPI_Comm grid_comm;

    dims[0] = Q;
    dims[1] = Q;
    periods[0] = 1;
    periods[1] = 1;
    reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm);
    MPI_Comm_rank(grid_comm, &grid_rank);
    MPI_Cart_coords(grid_comm, grid_rank, 2, coords);

    std::cout << grid_rank << " ; " << coords[0] << ':' << coords[1] << std::endl;

    // define sub-matrix type (?)

    int* submatrix = new int[N_over_Q*N_over_Q];
    MPI_Scatter(completeMatrix, N_over_Q*N_over_Q, MPI_INT, submatrix, N_over_Q*N_over_Q, MPI_INT, 0, grid_comm);

    std::string data = std::to_string(coords[0]) + ":" + std::to_string(coords[1]) + ":\t";
    for (int y = 0; y < N_over_Q; ++y) {
        for (int x = 0; x < N_over_Q; ++x) {
            data += std::to_string(submatrix[y*N_over_Q + x]);
            data += " ";
        }
        data += " ; ";
    }
    std::cout << data << std::endl;


    // run fox's algorithm until D > N on all nodes
        // for q steps
            // select starting nodes on diagonal, move diagonal by #step to the right
            // broadcast matrix to complete row, multiply received matrices with own data
            // send matrix to node directly above, multiply received matrix with own data

            // NOTE: use min-plus-matrix multiplication: choose minimum > 0 after element-wise addition


    // gather sub-matrices into matrix
    // print result / compare to to optional output file

    MPI_Comm_free(&grid_comm);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
