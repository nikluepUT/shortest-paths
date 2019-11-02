#include <iostream>

constexpr int INF = 999999;

void min_plus_matrix_multiply(int* own, const int* other, const int N){
    const auto MATRIX_SIZE = N * N;
    int* result = new int[MATRIX_SIZE];
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        result[i] = INF;
    }


    int min = 0;
    int tmp = 0;
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        min = own[i] + other[i];

        for (int k = 0; k < N; ++k) {
            tmp = own[(i/N)*N + k] + other[k*N + (i % N)];
            if (tmp < min){
                min = tmp;
            }
        }
        result[i] = min;
    }


    std::copy(result, result + MATRIX_SIZE, own);
    // cleanup
    delete[] result;
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
        min_plus_matrix_multiply(completeMatrix, completeMatrix, N);
    }

    print_matrix(completeMatrix, N);

    delete[] completeMatrix;
    return 0;
}