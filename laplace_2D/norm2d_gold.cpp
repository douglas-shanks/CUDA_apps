#include <math.h>       /* sqrt */
#include <stdio.h>      /* printf */
#include <iostream>
////////////////////////////////////////////////////////////////////////////////
double Gold_norm2d(int NX, int NY, double* u1, double* u2) 
{
  int   i, j, ind;
  double norm;
  
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {   // i loop innermost for sequential memory access
        ind = i + j*NX;
        norm += (u1[ind]-u2[ind])*(u1[ind]-u2[ind]);
      }
    }
        norm = sqrt(norm/ (double)(NX*NY));
        return norm;
}
