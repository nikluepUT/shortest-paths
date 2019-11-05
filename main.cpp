#include <iostream>
#include <cmath>
#include <algorithm>

#include <mpi.h>

constexpr int MASTER = 0;
constexpr int INF = std::numeric_limits<int>::max() / 2;

void min_plus_matrix_multiply(const int* matrixA, const int* matrixB, int* result, const int N);
void print_matrix(const int* matrix, const int N);

int main(int argc, char *argv[]) {

    // setup
    double startTime, endTime;
    int N, P, Q, N_OVER_Q, MATRIX_SIZE;
    int *completeMatrix = nullptr, *myMatrix, *result, *matrixA, *matrixB, *otherMatrixA;

    MPI_Comm commGrid, commRow;
    int myWorldRank, myGridRank, coords[2], upRank, downRank;

    std::cout << std::fixed;
    std::cout.precision(2);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &myWorldRank);

    // read input data, perform validation checks
    if (myWorldRank == MASTER) {

        std::cin >> N;

        // number of processes (P) must be a square, size of myMatrix (N) must be divisible by Q
        Q = static_cast<int>( sqrt(P) );
        if ((Q * Q) != P || N % Q != 0) {
            std::cout << "Incompatible values for P and N!" << std::endl;
            return EXIT_FAILURE;
        }

        // broadcast dimension data to other processes
        N_OVER_Q = N / Q;
        MATRIX_SIZE = N_OVER_Q*N_OVER_Q;
        MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        MPI_Bcast(&Q, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

        // read myMatrix in structurally more useful manner
        completeMatrix = new int[N * N];
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                int i, j, k;
                i = y/N_OVER_Q;
                j = x/N_OVER_Q;
                k = (y % N_OVER_Q)*N_OVER_Q + (x % N_OVER_Q);

                int input;
                std::cin >> input;
                if (input == 0 && x != y) {
                    input = INF;
                }
                completeMatrix[
                        i*N_OVER_Q*N +
                        j*MATRIX_SIZE +
                        k] = input;
            }
        }

        // start timing
        MPI_Barrier(MPI_COMM_WORLD);
        startTime = MPI_Wtime();

    } else {
        // receive dimension data
        MPI_Bcast(&N, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        MPI_Bcast(&Q, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        N_OVER_Q = N / Q;

        MPI_Barrier(MPI_COMM_WORLD);
    }
    MATRIX_SIZE = N_OVER_Q*N_OVER_Q;


    // create cartesian grid communicator, split into rows and gather info about own position / neighbours
    const int dims[] = {Q, Q}, periods[] = {1, 1}, includedDims[] = {0, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &commGrid);
    MPI_Cart_sub(commGrid, includedDims, &commRow);


    MPI_Comm_rank(commGrid, &myGridRank);
    MPI_Cart_coords(commGrid, myGridRank, 2, coords);

    coords[0] -= 1;
    MPI_Cart_rank(commGrid, coords, &upRank);
    coords[0] += 2;
    MPI_Cart_rank(commGrid, coords, &downRank);
    coords[0] -= 1;


    // scatter matrix

    result = new int[MATRIX_SIZE];
    myMatrix = new int[MATRIX_SIZE];
    matrixA = new int[MATRIX_SIZE];
    matrixB = new int[MATRIX_SIZE];
    otherMatrixA = new int[MATRIX_SIZE];

    MPI_Scatter(
            completeMatrix, MATRIX_SIZE, MPI_INT,
            myMatrix, MATRIX_SIZE, MPI_INT,
            MASTER, commGrid
    );

    std::copy(myMatrix, myMatrix + MATRIX_SIZE, matrixA);
    std::copy(myMatrix, myMatrix + MATRIX_SIZE, matrixB);
    std::fill(result, result + MATRIX_SIZE, INF);

    // run fox's algorithm until D > N on all nodes
    for (int D = 1; D < N; D *= 2) {

        // for q steps
        for (int step = 0; step < Q; ++step) {

            // select starting nodes on diagonal, move diagonal by #step to the right
            // broadcast myMatrix to complete row
            auto activeRank = (step + coords[0]) % Q;
            auto isActive = coords[1] == activeRank;

            if (isActive) {
                MPI_Bcast(matrixA, MATRIX_SIZE, MPI_INT, activeRank, commRow);
                min_plus_matrix_multiply(matrixA, matrixB, result, N_OVER_Q);
            } else {
                MPI_Bcast(otherMatrixA, MATRIX_SIZE, MPI_INT, activeRank, commRow);
                min_plus_matrix_multiply(otherMatrixA, matrixB, result, N_OVER_Q);
            }

            // send myMatrix to node directly above, multiply received myMatrix with own data
            MPI_Sendrecv_replace(matrixB, MATRIX_SIZE, MPI_INT, upRank, 0, downRank, 0, commGrid, nullptr);
        }

        std::copy(result, result + MATRIX_SIZE, matrixA);
        std::copy(result, result + MATRIX_SIZE, matrixB);
        std::fill(result, result + MATRIX_SIZE, INF);
    }


    // gather sub-matrices into myMatrix
    MPI_Gather(matrixA, MATRIX_SIZE, MPI_INT, completeMatrix, MATRIX_SIZE, MPI_INT, MASTER, MPI_COMM_WORLD);

    // print result / compare to to optional output file, display computation time
    if (myWorldRank == MASTER) {

        endTime = MPI_Wtime();
        std::cout << "Computations complete after " << endTime - startTime << " s:" << std::endl;

        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                int i, j, k;
                i = y / N_OVER_Q;
                j = x / N_OVER_Q;
                k = (y % N_OVER_Q) * N_OVER_Q + (x % N_OVER_Q);

                auto num = completeMatrix[
                        i * N_OVER_Q * N +
                        j * MATRIX_SIZE +
                        k];
                std::cout << (num >= INF ? 0 : num);
                std::cout << " ";
            }
            std::cout << std::endl;
        }
    }
    // cleanup
    if (myWorldRank == MASTER) {
        delete[] completeMatrix;
    }
    delete[] myMatrix;
    delete[] result;
    delete[] matrixA;
    delete[] matrixB;
    delete[] otherMatrixA;

    MPI_Comm_free(&commGrid);
    MPI_Comm_free(&commRow); // not sure if necessary, example does not include this
    MPI_Finalize();


    return EXIT_SUCCESS;
}

/**
 * The min-plus product C of two matrices A and B is defined as:
 * c_ij = min_k{a_ik + b_kj}
 * In other words: c_ij is the smallest sum of A k's row and B k's column.
 * The resulting product is copied to the own matrix.
 *
 * @param matrixA vector representing a matrix with dimension nxn
 * @param matrixB vector representing a matrix with dimension nxn
 * @param N matrix dimension
 */
void min_plus_matrix_multiply(const int* matrixA, const int* matrixB, int* result, const int N){
    const auto MATRIX_SIZE = N * N;


    int tmp = 0;
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int k = 0; k < N; ++k) {
            tmp = matrixA[(i / N) * N + k] + matrixB[k * N + (i % N)];
            if (tmp < result[i]){
                result[i] = tmp;
            }
        }
    }
}

void print_matrix(const int* matrix, const int N){
    for (int i = 0; i < N * N; ++i) {
        if (i > 0  && i % N == 0){
            std::cout << "\n";
        }
        std::cout << matrix[i] << " ";
    }
    std::cout << std::endl;
}
