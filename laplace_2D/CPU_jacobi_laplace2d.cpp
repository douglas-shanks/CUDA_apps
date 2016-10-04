
////////////////////////////////////////////////////////////////////////////////
void CPU_jacobi_laplace2d(int NX, int NY, double* u1, double* u2, double* b, double* alpha) 
{
  int   i, j, ind;
  double quart=1.0/4.0;  // predefining this improves performance more than 10%
  double hx, hy, hx2, hy2;
  
  hx = 1 / double (NX-1);
  hy = 1 / double (NY-1);
  hx2 = hx*hx;
  hy2 = hy*hy;
  
    for (j=0; j<NY; j++) {
      for (i=0; i<NX; i++) {   // i loop innermost for sequential memory access
		ind = i + j*NX;

        if (i==0 || i==NX-1) {
          u2[ind] = u1[ind]/hx2;          // Dirichlet b.c.'s
        }
        else if (j==0 || j==NY-1){
          u2[ind] = u1[ind]/hy2;          // Dirichlet b.c.'s	
        }
        else {
        
        u2[ind] =( u1[ind-1    ]/hx2 + u1[ind+1    ]/hx2 + u1[ind-NX   ]/hy2 +  u1[ind+NX   ]/hy2  + b[ind] )/( 2.0/hx2 + 2.0/hy2 ); 
        
         // u2[ind] = ( ( (alpha[ind-1] + alpha[ind])/2 )*u1[ind-1    ]/hx2 + ( (alpha[ind+1] + alpha[ind])/2 )*u1[ind+1    ]/hx2
         //           + ( (alpha[ind-NX] + alpha[ind - NX - 1])/2 )*u1[ind-NX   ]/hy2 +  ( (alpha[ind + NX] + alpha[ind + NX -1 ])/2 )*u1[ind+NX   ]/hy2  + b[ind]) /   ( (1.0 + alpha[ind])/hy2 + (1 + alpha[ind])/hx2 ); 
        }
      }
    }

}

