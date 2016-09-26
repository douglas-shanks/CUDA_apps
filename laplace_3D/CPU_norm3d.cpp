#include <math.h>       /* sqrt */
#include <stdio.h>      /* printf */
#include <iostream>
////////////////////////////////////////////////////////////////////////////////
double CPU_norm3d(int NX, int NY, int NZ, double* u1, double* u2) 
{
  int   i, j, k, ind;
  double norm;
  
  for (k=0; k<NZ; k++) {
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {   // i loop innermost for sequential memory access
        ind = i + j*NX+ k*NX*NY;
        norm += (u1[ind]-u2[ind])*(u1[ind]-u2[ind]);
      }
    }
  }  
        norm = sqrt(norm/ (double)(NX*NY*NZ));
        return norm;
}
