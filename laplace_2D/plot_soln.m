function dummy = plot_soln(m,inp)

fid = fopen(inp);

u = fscanf(fid,'%g');

fclose(fid);
h = 1/(m-1);
% we don't bother with the Dirichlet bndry
x = (h*[0:m-1]');
u = reshape(u,m,m);
surf(x,x,u);shading interp;
