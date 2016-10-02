% Grid points
NX = 40;
% max iterations
maxit = 500;
% output tolerance
tol = 1e-3;
% rhs, b = 1;
b = zeros(NX,NX);
% Solution
u = zeros(NX,NX);
% initial guess
uI = zeros(NX,NX);

h = 1/(NX-1);

for i = 1:NX
  for j = 1:NX
    if (i==1 || i==NX || j==1 || j==NX)
      uI(i,j) = 0.0;
    else
      uI(i,j) = 0.01;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end  

for i = 1:NX
  for j = 1:NX
    ind = i+ (j-1)*(NX);
    if (i==1 || i==NX || j==1 || j==NX)
      b(i,j) = 0.0;
    else
      b(i,j) = 1.0;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end  

% Iteration

for (k = 1: maxit)

  for (i = 1:NX)
    for (j = 1:NX)
      if (i==1 || i==NX || j==1 || j==NX)
        u(i,j) = uI(i,j);
      else   

        u(i,j) = ( uI(i-1,j) + uI(i+1,j) + uI(i,j-1) + uI(i,j+1) + h*h*b(i,j)  )*(1/4);
      end
    end
   end  
   
 temp = uI;
 uI = u;
 u = temp; 
      
end

