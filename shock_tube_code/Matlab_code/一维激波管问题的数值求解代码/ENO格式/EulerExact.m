function [rho,ux,p]=EulerExact(rho1,u1,p1,rho4,u4,p4,x,tEnd)
% function [rho,ux,p,e,t,Mach,entro]=EulerExact(rho1,u1,p1,rho4,u4,p4,x,tEnd)
%解析计算Riemann问题
%程序来源网络，用牛顿法求解解析解
%由于已经确定是Riemann问题，所以在计算的时候只需要给定左右初值

global gamma;
alpha=(gamma+1)/(gamma-1);

PRL = p4/p1;
cright = sqrt(gamma*p4/rho4); 
cleft  = sqrt(gamma*p1/rho1);
CRL = cright/cleft;
MACHLEFT = (u1-u4)/cleft;

% Basic shock tube relation equation (10.51)
f = @(P) (1+MACHLEFT*(gamma-1)/2-(gamma-1)*CRL*(P-1)/sqrt(2*gamma*(gamma-1+(gamma+1)*P)))^(2*gamma/(gamma-1))/P-PRL;

% solve for P = p34 = p3/p4
p34 = fzero(f,3);

p3 = p34*p4;
rho3 = rho4*(1+alpha*p34)/(alpha+p34); 
rho2 = rho1*(p34*p4/p1)^(1/gamma);
u2 = u1-u4+(2/(gamma-1))*cleft*(1-(p34*p4/p1)^((gamma-1)/(2*gamma)));
c2 = sqrt(gamma*p3/rho2);
spos = 0.5+tEnd*cright*sqrt((gamma-1)/(2*gamma)+(gamma+1)/(2*gamma)*p34)+tEnd*u4;

x0 = 0.5*(x(1)+x(end));
conpos=x0 + u2*tEnd+tEnd*u4;	% Position of contact discontinuity
pos1 = x0 + (u1-cleft)*tEnd;	% Start of expansion fan
pos2 = x0 + (u2+u4-c2)*tEnd;	% End of expansion fan

% Plot structures
p = zeros(size(x)); 
ux= zeros(size(x)); 
rho = zeros(size(x));
Mach = zeros(size(x));  
cexact = zeros(size(x));

for i = 1:length(x)
    if x(i) <= pos1
        p(i) = p1;
        rho(i) = rho1;
        ux(i) = u1;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = ux(i)/cexact(i);
    elseif x(i) <= pos2
        p(i) = p1*(1+(pos1-x(i))/(cleft*alpha*tEnd))^(2*gamma/(gamma-1));
        rho(i) = rho1*(1+(pos1-x(i))/(cleft*alpha*tEnd))^(2/(gamma-1));
        ux(i) = u1 + (2/(gamma+1))*(x(i)-pos1)/tEnd;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = ux(i)/cexact(i);
    elseif x(i) <= conpos
        p(i) = p3;
        rho(i) = rho2;
        ux(i) = u2+u4;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = ux(i)/cexact(i);
    elseif x(i) <= spos
        p(i) = p3;
        rho(i) = rho3;
        ux(i) = u2+u4;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = ux(i)/cexact(i);
    else
        p(i) = p4;
        rho(i) = rho4;
        ux(i) = u4;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = ux(i)/cexact(i);
    end
end
entro = log(p./rho.^gamma);	% entropy
e = p./((gamma-1).*rho);	% internal energy
t = 2/(2/(gamma-1)).*e;                 % temperature
end