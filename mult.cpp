#include <iostream>
#include <limits>


constexpr int INF = std::numeric_limits<int>::max() / 2;

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



int main() {
    int N;
    std::cin >> N;

    auto completeMatrix = new int[N * N];
    auto result = new int[N*N];
    std::fill(result, result + N*N, INF);
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {

            int input;
            std::cin >> input;
            if (input == 0 && x != y) {
                input = INF;
            }
            completeMatrix[y*N + x] = input;
        }
    }

    for (int D = 1; D < N; D *= 2) {
        min_plus_matrix_multiply(completeMatrix, completeMatrix, result, N);
        std::copy(result, result + N*N, completeMatrix);
    }

    print_matrix(result, N);

    delete[] completeMatrix;
    delete[] result;
    return 0;
}