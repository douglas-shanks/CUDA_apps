% Grid points
NX = 40;
% max iterations
maxit = 150;
% output tolerance
tol = 1e-3;
% rhs, b = 1;
b = zeros(NX*NX,1);
% Solution
u = zeros(NX*NX,1);
% initial guess
uI = zeros(NX*NX,1);
% alpha
alpha = ones(NX*NX,1);

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


% diffusion coefficient

for i = 1:NX
  for j = 1:NX
    ind = i+ (j-1)*(NX);
    if (i>=1 && i<=NX/2 && j>=1 || j<=NX/2)
      alpha(ind,1) = 200.0;
    else
      alpha(ind,1) = 1.0;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end 

% Take the harmonic mean
alphax = ones(NX*NX,1);
alphay = ones(NX*NX,1);
for i = 1:NX
  for j = 1:NX
      ind = i+ (j-1)*(NX);
      if (i==1 || i==NX || j==1 || j==NX )
      alphax(ind,1) = alpha(ind,1);
      alphay(ind,1) = alpha(ind,1);
      else
      alphax(ind,1) = 2*(alpha(ind,1)*alpha(ind+1,1))/(alpha(ind,1)+alpha(ind+1,1));
      alphay(ind,1) = 2*(alpha(ind,1)*alpha(ind+NX,1))/(alpha(ind,1)+alpha(ind+NX,1));
      end
  end
end 

% Iteration

for (k = 1: maxit)

  for (i = 1:NX)
    for (j = 1:NX)
    ind = i + (j-1)*NX;    
      if (i==1 || i==NX || j==1 || j==NX)
        u(ind) = uI(ind)/h/h;
      else
        u(ind) = ( alphax(ind,1)*uI(ind-1,1)/h/h + alphax(ind+1,1)*uI(ind+1,1)/h/h + alphay(ind)*uI(ind - NX,1)/h/h + alphay(ind+NX,1)*uI(ind + NX,1)/h/h + b(ind)  )/((alphax(ind,1) + alphax(ind+1,1))/h/h + (alphay(ind,1) + alphay(ind+NX,1))/h/h);
      end
    end
   end  
   
 temp = uI;
 uI = u;
 u = temp; 
      
end

