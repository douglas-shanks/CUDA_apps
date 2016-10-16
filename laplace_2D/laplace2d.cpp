//
// Program to solve Laplace equation on a regular 2D grid
// Just CPU code for now
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>

#define PI 4.0*atan(1.0)
////////////////////////////////////////////////////////////////////////
// declare Gold routines
////////////////////////////////////////////////////////////////////////

void CPU_jacobi_laplace2d(int NX, int NY, double* h_u1, double* h_u2, double* h_b, double* alpha);

double CPU_norm2d(int NX, int NY, double* h_u1, double* h_u2);

////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int main(int argc, const char **argv){


  // For timer
  std::clock_t start;
  double duration;

  // 'h_' prefix - CPU (host) memory space

  int    NX=40, NY=40, max_iters=5000, 
         bx, by, i, j, ind, iters=0;
  double  *h_u1, *h_u2, *h_u3, *h_b, *h_foo, *alpha, err, norm, tol = 1e-4, rx;

  printf("\nGrid dimensions: %d x %d\n", NX, NY);

  // allocate memory for arrays

  h_u1 = (double *)malloc(sizeof(double)*NX*NY);
  h_u2 = (double *)malloc(sizeof(double)*NX*NY);
  h_u3 = (double *)malloc(sizeof(double)*NX*NY);
  h_b = (double *)malloc(sizeof(double)*NX*NY);
  h_foo = (double *)malloc(sizeof(double)*NX*NY);
  alpha = (double *)malloc(sizeof(double)*NX*NY);

  // initialise u1. Initial guess for iteration, say u=1 which result in a bowl

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX;
        if (i==0 || i==NX-1 || j==0 || j==NY-1)
          h_u1[ind] = 0.0;           // Dirichlet b.c.'s
        else
          h_u1[ind] = 0.01;
      }
    }

  // initialise rhs b. Start with b = 0 

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX;
        if (i==0 || i==NX-1 || j==0 || j==NY-1)
          h_b[ind] = 0.0;           // Dirichlet b.c.'s
        else
          h_b[ind] = 1.0;//sin(PI*i/NX)*sin(PI*j/NY);
      }
    }
  // initialise alpha in D(alpha D u) = 0 

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
          ind = i + j*NX;
        
        //if( (round(0.5*NX)<=i && i<=round(0.75*NX)) && (0.0<=j && j<=round(0.25*NY)))
   		//	alpha[ind]=10.0;
		//else if((round(0.25*NX)<=i && i<=round(0.5*NX))&&(round(0.25*NY)<=j && j<=round(0.5*NY)))
   		//	alpha[ind]=10.0;
		//else if((round(0.5*NX)<=i && i<=round(0.75*NX))&&(round(0.5*NY)<=j && j<=round(0.75*NY)))
   		//	alpha[ind]=10.0;
		//else if((round(0.25*NX)<=i && i<=round(0.5*NX))&&(round(0.75*NY)<=j && j<=round(1.0*NY)))
   		//	alpha[ind]=10.0;
		//else
   			//alpha[ind] = 0.001;
        
          alpha[ind] = 1.0;
      }
    }
    
 // print the solution to a text file
std::ofstream out("alpha.txt");
for (j=0; j<NY; j++) {
	for (i=0; i<NX; i++) {
        ind = i + j*NX;
        out << alpha[ind]<< '\n';
    }
}
out.close();   
    
  // Jacobi iteraton

  start = std::clock();
  norm = 0.0;
  for (int i = 1; i <= max_iters; ++i) {
  
    CPU_jacobi_laplace2d(NX, NY, h_u1, h_u3, h_b, alpha);
    
    // Add a declaration to exit the loop if the tolerance is reached.
    norm = CPU_norm2d(NX, NY, h_u1, h_u3);  
     
    iters = iters + 1; //iteration counter
    if (norm <= tol) break;
    
    h_foo = h_u1; h_u1 = h_u3; h_u3 = h_foo;   // swap h_u1 and h_u3
  }
  
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  printf("\n%d iterations of Jacobi: %.8f (s) \n \n", iters, duration);
  printf("2-norm of residual= %.16f \n",norm);

  // error check

  err = 0.0;

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {
        ind = i + j*NX;
        err += (h_u1[ind]-h_u2[ind])*(h_u1[ind]-h_u2[ind]);
      }
    }

  printf("rms error = %.16f \n",sqrt(err/ (double)(NX*NY)));

// print the solution to a text file
std::ofstream outlap("laplace.txt");
for (j=0; j<NY; j++) {
	for (i=0; i<NX; i++) {
        ind = i + j*NX;
        outlap << h_u1[ind]<< '\n';
    }
}
out.close();

 // Release CPU memory

  free(h_u1);
  free(h_u2);
  free(h_u3);
  free(h_b);

}
