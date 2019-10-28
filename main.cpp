#include <iostream>

int main() {

    // read N and P, compute Q, validate that all values are compatible with Fox's algorithm

    // read matrix

    // broadcast N/Q
    // create cartesian grid communicator
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

    return 0;
}