function dF = FluxDiffLax_W(U,F,dx,dt)
global gamma;

%基础差分的部分
dF = 0.5*(F(:,3:end) - F(:,1:end-2))/dx;

%二阶导数的修正
%计算Jocabi矩阵
AR = zeros(3,3,size(F,2)-2);
UR = 0.5*(U(:,2:end-1)+U(:,3:end));
AR(1,1,:)=0;  
AR(1,2,:)=1;  
AR(1,3,:)=0;
AR(2,1,:) = -(3-gamma)/2 * UR(2,:).^2./UR(1,:).^2;
AR(2,2,:) = (3-gamma)*UR(2,:)./UR(1,:);
AR(2,3,:) = gamma-1;
AR(3,1,:) = -gamma*UR(2,:).*UR(3,:)./UR(1,:).^2 + (gamma-1)*UR(2,:).^3./UR(1,:).^3;
AR(3,2,:) = gamma*UR(3,:)./UR(1,:) - 3/2*(gamma-1)*UR(2,:).^2./UR(1,:).^2;
AR(3,3,:) = gamma*UR(2,:)./UR(1,:);

AL = zeros(3,3,size(F,2)-2);
UL = 0.5*(U(:,1:end-2)+U(:,2:end-1));
AL(1,1,:)=0;  
AL(1,2,:)=1;  
AL(1,3,:)=0;
AL(2,1,:) = -(3-gamma)/2 * UL(2,:).^2./UL(1,:).^2;
AL(2,2,:) = (3-gamma)*UL(2,:)./UL(1,:);
AL(2,3,:) = gamma-1;
AL(3,1,:) = -gamma*UL(2,:).*UL(3,:)./UL(1,:).^2 + (gamma-1)*UL(2,:).^3./UL(1,:).^3;
AL(3,2,:) = gamma*UL(3,:)./UL(1,:) - 3/2*(gamma-1)*UL(2,:).^2./UL(1,:).^2;
AL(3,3,:) = gamma*UL(2,:)./UL(1,:);
%计算修正,为了方便阅读没有进行向量化
dF2 = zeros(3,size(F,2)-2);
for i = 2:size(dF2,2)+1
    ARi = AR(:,:,i-1);
    ALi = AL(:,:,i-1);
    mid_R = 0.5*dt/dx^2 * ARi*(F(:,i+1)-F(:,i));
    mid_L = 0.5*dt/dx^2 * ALi*(F(:,i)-F(:,i-1));
    dF2(:,i-1) = mid_R - mid_L;
end

%最终结果
dF = dF - dF2;


end