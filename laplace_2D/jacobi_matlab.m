% Grid points
NX = 40;
% max iterations
maxit = 500;
% output tolerance
tol = 1e-3;
% rhs, b = 1;
b = zeros(NX*NX,1);
% Solution
u = zeros(NX*NX,1);
% initial guess
uI = zeros(NX*NX,1);

h = 1/(NX-1);

for i = 1:NX
  for j = 1:NX
    ind = i+ (j-1)*(NX);
    if (i==1 || i==NX || j==1 || j==NX)
      uI(ind,1) = 0.0;
    else
      uI(ind,1) = 0.01;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end  

for i = 1:NX
  for j = 1:NX
    ind = i+ (j-1)*(NX);
    if (i==1 || i==NX || j==1 || j==NX)
      b(ind,1) = 0.0;
    else
      b(ind,1) = 1.0;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end  

% Iteration

for (k = 1: maxit)

  for (i = 1:NX)
    for (j = 1:NX)
    ind = i + (j-1)*NX;    
      if (i==1 || i==NX || j==1 || j==NX)
        u(ind) = uI(ind);
      else
        u(ind) = ( uI(ind-1) + uI(ind+1) + uI(ind - NX) + uI(ind + NX) + h*h*b(ind)  )*(1/4);
      end
    end
   end  
   
 temp = uI;
 uI = u;
 u = temp; 
      
end

