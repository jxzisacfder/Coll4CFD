function dF = SpaceDiscre(U,dt,dx,x)
%Roe格式的空间分裂
global gamma;

%首先计算通量，Roe格式不需要进行通量分裂
[rho,u,p] = conservation_to_physics(U);
E = p/(gamma-1) + 0.5*rho.*u.^2;
%特征值
if ( sum(p./rho <= 0) > 0 )
    [p;rho]
    error('出现负数');
end
c = sqrt(gamma*p./rho);
lambda = max(max(abs(u)+c));
%通量
% F = [rho.*u; rho.*u.^2+p; (E+p).*u];
% F_posi = 0.5*(F+lambda*U);
% F_nega = 0.5*(F-lambda*U);

%计算正负通量的迎风差分，Roe格式在内部计算通量F，但是尽量保持一致的程序框架
dF = FluxDiffRoe(U,dx,dt);

end
