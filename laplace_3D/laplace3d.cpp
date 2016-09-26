//
// Program to solve Laplace equation on a regular 3d grid
// Just CPU code for now
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <cstdio>
#include <ctime>

////////////////////////////////////////////////////////////////////////
// declare Gold routine
////////////////////////////////////////////////////////////////////////

void CPU_Jacobi3d(int NX, int NY, int NZ, double* h_u1, double* h_u2, double* h_b);

double CPU_norm3d(int NX, int NY, int NZ, double* h_u1, double* h_u2);

////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv){


  // For timer
  std::clock_t start;
  double duration;

  // 'h_' prefix - CPU (host) memory space

  int    NX=150, NY=150, NZ=150, max_iters=1000, 
         bx, by, i, j, k, ind, iters=0;
  double  *h_u1, *h_u2, *h_u3, *h_b, *h_foo, err, norm, tol = 1e-2;

  printf("\nGrid dimensions: %d x %d x %d\n", NX, NY, NZ);

  // allocate memory for arrays

  h_u1 = (double *)malloc(sizeof(double)*NX*NY*NZ);
  h_u2 = (double *)malloc(sizeof(double)*NX*NY*NZ);
  h_u3 = (double *)malloc(sizeof(double)*NX*NY*NZ);
  h_b = (double *)malloc(sizeof(double)*NX*NY*NZ);

  // initialise u1

  for (k=0; k<NZ; k++) {
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX + k*NX*NY;

        if (i==0 || i==NX-1 || j==0 || j==NY-1|| k==0 || k==NZ-1)
          h_u1[ind] = 1.0;           // Dirichlet b.c.'s
        else
          h_u1[ind] = 0.0;
      }
    }
  }
  
  // initialise rhs b. Start with b = 0 and solve the error equation
  
  for (k=0; k<NZ; k++) {
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX + k*NX*NY;

        if (i==0 || i==NX-1 || j==0 || j==NY-1|| k==0 || k==NZ-1)
          h_b[ind] = 0.0;           // Dirichlet b.c.'s
        else
          h_b[ind] = 0.0;
      }
    }
  }
  // Jacobi iteraton

  start = std::clock();
  norm = 0.0;
  for (int i = 1; i <= max_iters; ++i) {
  
    CPU_Jacobi3d(NX, NY, NZ, h_u1, h_u3, h_b);
    
    // Add a declaration to exit the loop if the tolerance is reached.
    norm = CPU_norm3d(NX, NY, NZ, h_u1, h_u3);
    // what is the norm
    printf("norm = %.16f \n",norm);
    
    iters = iters + 1; //iteration counter
    if (norm <= tol) break;
    
    h_foo = h_u1; h_u1 = h_u3; h_u3 = h_foo;   // swap h_u1 and h_u3
  }

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  printf("\n%d iterations of Jacobi: %.8f (s) \n \n", iters, duration);
  printf("2-norm of residual = %.16f \n",norm);

  // error check

  err = 0.0;

  for (k=0; k<NZ; k++) {
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX + k*NX*NY;
        err += (h_u1[ind]-h_u2[ind])*(h_u1[ind]-h_u2[ind]);
      }
    }
  }

  printf("rms error = %f \n",sqrt(err/ (double)(NX*NY*NZ)));

 // Release CPU memory

  free(h_u1);
  free(h_u2);
  free(h_u3);
  free(h_b);

}
