[rho0,u0,p0] = Euler_IC(x,1);
[rho_exact,u_exact,p_exact] = EulerExact(rho0(1),u0(1),p0(1),rho0(end),u0(end),p0(end),x,i*dt);
E_exact = p_exact./((gamma-1)*rho_exact)+0.5*u_exact.^2;
rho_error = norm(rho-rho_exact)*dx;
u_error = norm(vel-u_exact)*dx;
p_error = norm(pre-p_exact)*dx;

%delete('test.csv');
flag = exist('test.csv','file');
fileID = fopen('test.csv','a');
if(flag~=2)
    fprintf(fileID,'%12s,%12s,%12s,%12s,%12s,%12s,%12s,%12s\n',...
        '#$\Delta x$','#$\Delta t$','#Error\_$\rho$','#Order\_$\rho$','#Error\_$u$','#Order\_$u$','#Error\_$p$','#Order\_$p$');
end
fprintf(fileID,'%.2E,%.2E,%.2E,%.2f,%.2E,%.2f,%.2E,%.2f\n',dx,dt,rho_error,0,u_error,0,p_error,0);
fclose(fileID);