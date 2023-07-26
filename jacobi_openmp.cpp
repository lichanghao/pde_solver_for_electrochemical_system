// Iteratively solve Poisson equation by Jacobi method
// Author: Changhao Li (changhaoli1997@gmail.com)
// Date created: 6/12/2023
// Date last modified: 7/11/2023

#include <ctime>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <array>
#include <omp.h>

using namespace std;

int main(int argc, char* argv[]) {
    cout << "Please enter the number of procs" << endl;
    int nproc = 1;
    cin >> nproc;
    omp_set_num_threads(nproc);
    
    // parameters
    int Lx = 4000, Ly = 4000;
    int max_iter = 5000;
    double delta = 1;
    double thres = 0.001;
    double prev_error = 1e-8;

    // set up 2D zero matrix u, b
    double **u = new double*[Ly+2];
    double **b = new double*[Ly+2];
    double **u_new = new double*[Ly+2];
    for (int i = 0; i < Ly+2; i++) {
        u[i] = new double[Lx+2];
        b[i] = new double[Lx+2];
        u_new[i] = new double[Lx+2];
    }
    for (int i = 0; i < Ly+2; i++) {
        for (int j = 0; j < Lx+2; j++) {
            u[i][j] = (double)(rand()) / (double)(RAND_MAX);
            b[i][j] = 0;
            u_new[i][j] = 0;
        }
    }

    // iteration 
    // Solving the Laplace equation \laplacian{u} = 0
    // which boils down to solving the linear equation Au = 0
    for (int iter = 0; iter < max_iter; iter++) {

        // Boundary condition
        for (int iy = 0; iy < Ly+2; iy++) {
            u[iy][0] = 0;
            u[iy][Lx+1] = 1;
        }
        for (int ix = 0; ix < Lx+2; ix++) {
            u[0][ix] = u[1][ix];
            u[Ly+1][ix] = u[Ly][ix];
        }
        
        // main loop for voxel-wise numerical iterations
        #pragma omp parallel for 
        #pragma& defaults(shared) private(iy, ix)
        for (int iy = 1; iy < Ly+1; iy++) {
            for (int ix = 1; ix < Lx+1; ix++) {
                u_new[iy][ix] = 0.25*(u[iy+1][ix] + u[iy-1][ix] + u[iy][ix+1] + u[iy][ix-1]) 
                                - 0.25*delta*delta*b[iy][ix];
            }
        }

        // calculate relative error
        double error = 0;
        #pragma omp parallel for 
        #pragma& defaults(shared) private(iy, ix)
        for (int iy = 1; iy < Ly+1; iy++) {
            for (int ix = 1; ix < Lx+1; ix++) {
                error += abs(0.25*(u[iy+1][ix] + u[iy-1][ix] + u[iy][ix+1] + u[iy][ix-1] - 4*u[iy][ix]) 
                         - 0.25*delta*delta*b[iy][ix]);
            }
        }
        double relative_error = abs(error - prev_error) / prev_error;
        prev_error = error;

        // update
        for (int iy = 1; iy < Ly+1; iy++) {
            for (int ix = 1; ix < Lx+1; ix++) {
                u[iy][ix] = u_new[iy][ix];
            }
        }
        if (relative_error <= thres) {
            cout << "Converged, total iterations: " << iter << endl;
            break;
        }

        // output
        if (iter % 20 == 0)
        cout << "relative error: " << relative_error << endl;

    }

    string outputfile1 = "solution.raw";
    ofstream file3(outputfile1.c_str(), ios::out | ios::binary);
    if (!file3) { cout << "Error writing "+ outputfile1 + "\n"; return 0; }
    for (int iy = 0; iy < Ly+2; iy++) {
        file3.write((char*)u[iy], (Lx+2)*sizeof(double));
    }
    file3.close();

    return 0;
}
