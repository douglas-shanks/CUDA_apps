% Grid points
NX = 40;
% max iterations
maxit = 150;
% output tolerance
tol = 1e-3;
% rhs, b = 1;
b = zeros(NX,NX);
% Solution
u = zeros(NX,NX);
% initial guess
uI = zeros(NX,NX);
% alpha
alpha = ones(NX,NX);

h = 1/(NX-1);

% initial guess
for i = 1:NX
  for j = 1:NX
    if (i==1 || i==NX || j==1 || j==NX)
      uI(i,j) = 0.0;
    else
      uI(i,j) = 0.01;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end  

% right hand side
for i = 1:NX
  for j = 1:NX
    if (i==1 || i==NX || j==1 || j==NX)
      b(i,j) = 0.0;
    else
      b(i,j) = 1.0;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end  

% diffusion coefficient

for i = 1:NX
  for j = 1:NX
    if (i>=1 && i<=NX/2 && j>=1 || j<=NX/2)
      alpha(i,j) = 2.0;
    else
      alpha(i,j) = 1.0;%sin(i*pi/NX)*sin(j*pi/NX);
    end
  end
end 

% Take the harmonic mean

for i = 1:NX
  for j = 1:NX
      if (i==1 || i==NX )
      alphax(i,j) = alpha(i,j);
      elseif (j==1 || j==NX )
      alphay(i,j) = alpha(i,j);
      else
      alphax(i,j) = 2*(alpha(i,j)*alpha(i+1,j))/(alpha(i,j)+alpha(i+1,j));
      alphay(i,j) = 2*(alpha(i,j)*alpha(i,j+1))/(alpha(i,j)+alpha(i,j+1));
      end
  end
end 

% Jacobi Iteration

for (k = 1: maxit)

  for (i = 1:NX)
    for (j = 1:NX)
      if (i==1 || i==NX || j==1 || j==NX)
        u(i,j) = uI(i,j)/h/h;
      else   
        u(i,j) = ( alphax(i,j)*uI(i-1,j)/h/h + alphax(i+1,j)*uI(i+1,j)/h/h + alphay(i,j)*uI(i,j-1)/h/h + alphay(i,j+1)*uI(i,j+1)/h/h + b(i,j)  )/((alphax(i+1,j)+alphax(i,j))/h/h + (alphay(i,j+1)+alphay(i,j))/h/h);
      end
    end
   end  
   
 temp = uI;
 uI = u;
 u = temp; 
      
end

