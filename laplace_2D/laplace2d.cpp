//
// Program to solve Laplace equation on a regular 2D grid
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
// declare Gold routines
////////////////////////////////////////////////////////////////////////

void Gold_laplace2d(int NX, int NY, double* h_u1, double* h_u2, double* h_b);

double Gold_norm2d(int NX, int NY, double* h_u1, double* h_u2);

////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv){


  // For timer
  std::clock_t start;
  double duration;

  // 'h_' prefix - CPU (host) memory space

  int    NX=125, NY=125, max_iters=10000, 
         bx, by, i, j, ind, iters=0;
  double  *h_u1, *h_u2, *h_u3, *h_b, *h_foo, err, norm, tol = 1e-2;

  printf("\nGrid dimensions: %d x %d\n", NX, NY);

  // allocate memory for arrays

  h_u1 = (double *)malloc(sizeof(double)*NX*NY);
  h_u2 = (double *)malloc(sizeof(double)*NX*NY);
  h_u3 = (double *)malloc(sizeof(double)*NX*NY);
  h_b = (double *)malloc(sizeof(double)*NX*NY);

  // initialise u1. Initial guess for iteration

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX;

        if (i==0 || i==NX-1 || j==0 || j==NY-1)
          h_u1[ind] = 1.0;           // Dirichlet b.c.'s
        else
          h_u1[ind] = 0.0;
      }
    }
    
  // initialise rhs b. Start with b = 0 and solve the error equation

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX;

        if (i==0 || i==NX-1 || j==0 || j==NY-1)
          h_b[ind] = 0.0;           // Dirichlet b.c.'s
        else
          h_b[ind] = 0.0;
      }
    }

  // Jacobi iteraton

  start = std::clock();
  norm = 0.0;
  for (int i = 1; i <= max_iters; ++i) {
    Gold_laplace2d(NX, NY, h_u1, h_u3, h_b);
    
    // Add a declaration to exit the loop if the tolerance is reached.
    norm = Gold_norm2d(NX, NY, h_u1, h_u3);
    // what is the norm
    printf("norm = %.16f \n",norm);
    if (norm <= tol) break;
    h_foo = h_u1; h_u1 = h_u3; h_u3 = h_foo;   // swap h_u1 and h_u3
    iters = iters + 1; //iteration counter
  }
  
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  printf("\n%d iterations of Jacobi: %.8f (s) \n \n", iters, duration);
  printf("2-norm = %.16f \n",norm);
  
  // print out corner of array

  /*
  for (k=0; k<3; k++) {
    for (j=0; j<8; j++) {
      for (i=0; i<8; i++) {
        ind = i + j*NX + k*NX*NY;
        printf(" %5.2f ", h_u1[ind]);
      }
      printf("\n");
    }
    printf("\n");
  }
  */

  // error check

  err = 0.0;

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX;
        err += (h_u1[ind]-h_u2[ind])*(h_u1[ind]-h_u2[ind]);
      }
    }

  printf("rms error = %.16f \n",sqrt(err/ (double)(NX*NY)));

 // Release CPU memory

  free(h_u1);
  free(h_u2);
  free(h_u3);

}
