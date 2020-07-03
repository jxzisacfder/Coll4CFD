%求解一维激波管问题
%ENO格式，计算格式中的光滑监测子并未使用最原始的检测子
%last modified in 2019.5.25
clc; clear;
format long e;

%计算参数设置
global N;
N = 1601;
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
U_last = physics_to_conservation(rho0,u0,p0);
for i = 1:Tnumber
    
    %Roe格式是空间离散格式，时间离散我们使用三阶RK3格式
    U_last = BC(U_last,3);
    U1 = U_last(:,4:end-3) - dt*SpaceDiscre(U_last,dt,dx,x);
    U1 = BC(U1,3);
    U2 = 0.75*U_last(:,4:end-3) + 0.25*U1(:,4:end-3) - 0.25*dt*SpaceDiscre(U1,dt,dx,x);
    U2 = BC(U2,3);
    U = 1/3*U_last(:,4:end-3) + 2/3*U2(:,4:end-3) - 2/3*dt*SpaceDiscre(U2,dt,dx,x);

    %维护循环并作图
    U_last = U;
    ResultPlot;
    pause(0.01);
end
cal_order;
