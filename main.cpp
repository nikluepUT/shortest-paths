#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <errno.h>
#include <cstring>
#include <cmath>

#define MASTER 0

char* FILEPATH = "/home/moritz/CLionProjects/shortest-paths/matrix_examples/input6";
int N;
int P;
int Q;
int N_over_Q;

int read_dimension(char *filename);
void read_matrix(char *filename, int** matrix);

int main(int argc, char *argv[]) {
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MASTER){
        // read N and P, compute Q, validate that all values are compatible with Fox's algorithm
        if (read_dimension(FILEPATH) < 0) {
            printf("Abort program.\n");
            return EXIT_FAILURE;
        }
        // read matrix
        read_matrix(FILEPATH, nullptr);
    }
    // broadcast N/Q
    MPI_Bcast(&N_over_Q, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // create cartesian grid communicator
    MPI_Comm grid_comm;
    int dim[2], period[2], reorder;

    dim[0] = Q;
    dim[1] = Q;
    period[0] = 0;
    period[1] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, 1, &grid_comm);

    // define sub-matrix type (?)

    // scatter matrix into sub-matrices


    // run fox's algorithm until D > N on all nodes
        // for q steps
            // select starting nodes on diagonal, move diagonal by #step to the right
            // broadcast matrix to complete row, multiply received matrices with own data
            // send matrix to node directly above, multiply received matrix with own data

            // NOTE: use min-plus-matrix multiplication: choose minimum > 0 after element-wise addition


    // gather sub-matrices into matrix
    // print result / compare to to optional output file

    MPI_Finalize();
    return EXIT_SUCCESS;
}

int read_dimension(char *filename){
    FILE *fp;

    // try to open file
    fp = fopen(filename, "r");
    if (fp == NULL){
        perror("fopen()");
        return -1;
    }

    // read dimensions and save it to global variable
    char* dim_line;
    size_t len = 0;
    if (getline(&dim_line, &len, fp) != -1){
        N = atoi(dim_line);
    }
    printf("N: %d \n",N);
    printf("P: %d \n",P);
    int Q = sqrt(P);
    if ((Q * Q) == P && (N % Q == 0)){
        printf("Q: %d \n",Q);
        N_over_Q = N/Q;
        printf("N/Q: %d \n", N_over_Q);
    }
    else{
        printf("No Q exists with P = Q * Q and N mod Q = 0\n");
        return -1;
    }

    if (dim_line)  {
        free(dim_line);
    }
    fclose(fp);

    return 0;
}

void read_matrix(char *filename, int** matrix) {
    FILE *fp;

    // try to open file
    fp = fopen(filename, "r");
    if (fp == NULL){
        perror("fopen()");
        return;
    }

    matrix = (int**) malloc(N * sizeof(int*));
    for (int i = 0; i < N; i++){
        matrix[i] = (int*) malloc (N * sizeof(int));
    }

    // read line by line and fill matrix
    char* line = NULL;
    size_t len = 0;
    int i = 0;
    getline(&line, &len, fp);
    while(getline(&line, &len, fp) != -1){
        //printf("%s", line);
        int j = 0;
        int start = 0;
        for (int end = 0; end < strlen(line); end++){
            if (line[end] == ' ' || line[end] == '\n'){
                char num[end-start + 1];
                memcpy(num, &line[start], end-start);
                num[end-start] = '\0';
                //printf("%s", num);

                start = end;

                matrix[i][j] = atoi(num);
                j++;
            }
        }
        i++;
    }

    printf("Matrix: \n");
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }

    if (line)  {
        free(line);
    }
    fclose(fp);
}