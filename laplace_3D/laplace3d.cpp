//
// Program to solve Laplace equation on a regular 3D grid
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

void Gold_laplace3d(int NX, int NY, int NZ, float* h_u1, float* h_u2);

////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv){


  // For timer
  std::clock_t start;
  double duration;

  // 'h_' prefix - CPU (host) memory space

  int    NX=150, NY=150, NZ=150, REPEAT=100, // Repeat = num of iterations
         bx, by, i, j, k, ind;
  float  *h_u1, *h_u2, *h_u3, *h_foo, err, norm;

  printf("\nGrid dimensions: %d x %d x %d\n", NX, NY, NZ);

  // allocate memory for arrays

  h_u1 = (float *)malloc(sizeof(float)*NX*NY*NZ);
  h_u2 = (float *)malloc(sizeof(float)*NX*NY*NZ);
  h_u3 = (float *)malloc(sizeof(float)*NX*NY*NZ);

  // initialise u1

  for (k=0; k<NZ; k++) {
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX + k*NX*NY;

        if (i==0 || i==NX-1 || j==0 || j==NY-1|| k==0 || k==NZ-1)
          h_u1[ind] = 1.0f;           // Dirichlet b.c.'s
        else
          h_u1[ind] = 0.0f;
      }
    }
  }


  // Jacobi iteraton

  start = std::clock();
  for (int i = 1; i <= REPEAT; ++i) {
    Gold_laplace3d(NX, NY, NZ, h_u1, h_u3);
    
    // Add a declaration to exit the loop if the tolerance is reached.
    // This will require
    
    h_foo = h_u1; h_u1 = h_u3; h_u3 = h_foo;   // swap h_u1 and h_u3
  }

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  printf("\n%dx Gold_laplace3d: %.8f (s) \n \n", REPEAT, duration);

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

  for (k=0; k<NZ; k++) {
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX + k*NX*NY;
        err += (h_u1[ind]-h_u2[ind])*(h_u1[ind]-h_u2[ind]);
      }
    }
  }

  printf("rms error = %f \n",sqrt(err/ (float)(NX*NY*NZ)));

 // Release CPU memory

  free(h_u1);
  free(h_u2);
  free(h_u3);

}
