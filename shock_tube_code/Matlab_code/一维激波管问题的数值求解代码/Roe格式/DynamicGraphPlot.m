%动态图作图
[rho,vel,pre] = conservation_to_physics(U);

if (mod(i,10)==0) 
    %计算当前时间下的解析解
    [rho_exact,u_exact,p_exact] = ...
             EulerExact(rho0(1),u0(1),p0(1),rho0(end),u0(end),p0(end),x,i*dt);
    E_exact = p_exact./((gamma-1)*rho_exact)+0.5*u_exact.^2;

    %作图
    figure(2);
    subplot(2,2,1);
    plot(x,rho,'b',x,rho_exact,'r');
    title('密度');
    subplot(2,2,2);
    plot(x,vel,'b',x,u_exact,'r');
    title('速度');
    subplot(2,2,3);
    plot(x,pre,'b',x,p_exact,'r');
    title('压力'); 

    xlabel(strcat('time = ',num2str(dt*i)));
    pause(0.01);
end
