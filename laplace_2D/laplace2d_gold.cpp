
////////////////////////////////////////////////////////////////////////////////
void Gold_laplace2d(int NX, int NY, double* u1, double* u2, double* b) 
{
  int   i, j, ind;
  double quart=1.0/4.0;  // predefining this improves performance more than 10%

    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {   // i loop innermost for sequential memory access
	ind = i + j*NX;

        if (i==0 || i==NX-1 || j==0 || j==NY-1) {
          u2[ind] = u1[ind] + b[ind];          // Dirichlet b.c.'s
        }
        else {
          u2[ind] = ( u1[ind-1    ] + u1[ind+1    ]
                    + u1[ind-NX   ] + u1[ind+NX   ]  + b[ind]) * quart;
        }
      }
    }

}

