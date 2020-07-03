function U = BC(U,n)
%边界处理的函数
%U是守恒量，n代表左右各延拓的长度

% %周期边界
% Uleft = U(:,end-n+1:end);
% Uright = U(:,1:n);
% U = [Uleft U Uright];

%输运边界
for i = 1:n
    U = [U(:,1) U U(:,end)];
end 

% %反射边界
% [rho,u,p] = conservation_to_physics(U);
% %边界直接赋值的部分
% u(1) = 0;      u(end) = 0;
% p(1) = p(2);   p(end-1) = p(end);
% %虚拟点的部分
% uleft = []; uright = [];
% pleft = []; pright = [];
% rholeft = []; rhoright = [];
% for i = 1:n
%     uleft = [-u(i+1) uleft];
%     uright = [uright -u(end-i)];
%    
%     pleft = [p(i+1) pleft];
%     pright = [pright p(end-i)];
%     
%     rholeft = [rho(i+1) rholeft];
%     rhoright = [rhoright rho(end-i)];    
% end
% u = [uleft u uright];
% p = [pleft p pright];
% rho = [rholeft rho rhoright];
% U = physics_to_conservation(rho,u,p);

end