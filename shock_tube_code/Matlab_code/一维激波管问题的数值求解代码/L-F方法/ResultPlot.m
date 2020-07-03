if mod(i*dt,0.01)==0
    figure(1);
    [rho,vel,pre] = conservation_to_physics(U);
    p1 = plot(x,rho,'b',x,vel,'r',x,pre,'black');set(p1,'linewidth', 3);
    grid on;
    legend('Density','Velocity','Pressure');
    title(['Lax\_Friedrichs scheme:  t = ' num2str(i*dt)]);
    drawnow;
end