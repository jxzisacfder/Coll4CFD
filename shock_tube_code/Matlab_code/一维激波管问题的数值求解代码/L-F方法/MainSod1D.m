%求解一维激波管问题
%last modified in 2019.5.24
clc; clear;
format long e;

%计算参数设置
N = 1001;
T = 0.10;
x = linspace(0,1,N);
dx = x(2) - x(1);
dt = dx*0.1;
Tnumber = T/dt;

%物理参数
global gamma;
gamma = 1.4;

%设置IC
[rho0,u0,p0] = Euler_IC(x,1);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % 总能量密度
a0 = sqrt(gamma*p0./rho0);            % 声速

%开始计算
%由于L-F格式为中心对称格式，不需要进行通量分裂
U_last = physics_to_conservation(rho0,u0,p0);
for i = 1:Tnumber
    
    %先处理边界，然后做空间离散，最后随时间推进
    U_last = BC(U_last,1); %当前问题使用输运边界即可
    dF = SpaceDiscre(U_last,dt,dx,x);
    U = 0.5*(U_last(:,1:end-2)+U_last(:,3:end))-dt*dF;
    
    %维护循环并作图
    U_last = U;
    DynamicGraphPlot;
%     pause(0.01);
    
end
