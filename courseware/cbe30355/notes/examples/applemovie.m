%% Movie of apple drainage
% We create a movie describing the drainage of a thin film around a sphere.
% The problem has been rendered dimensionless, so we start with an
% initially uniform film depth of unity and evaluate it as a function of
% theta.

for k=1:50 %repeat the movie 50 times
n=100;
theta=[0:n]/n*pi;
delta=ones(size(theta));

dt=0.0001;

t=0;

while t<0.25
    deltadot=zeros(size(delta));
    for i=2:n
        deltadot(i)=-2*cos(theta(i))*delta(i)^3-3*sin(theta(i))*delta(i)^2*(delta(i+1)-delta(i-1))*pi/n;
    end
    deltadot(1)=-2*delta(1)^3;
    deltadot(n+1)=2*delta(n+1)^3;
    
    delta=delta+deltadot*dt;
    
    t=t+dt;
    
    figure(1)
    subplot(1,2,1), plot(theta,delta)
    xlabel('theta')
    ylabel('delta')
    title('Sphere Drainage Thickness')
%    title(['Dimensionless thickness at t = ',num2str(t)])
    
    x0=sin(theta);
    y0=cos(theta);
    x=(1+.1*delta).*sin(theta);
    y=(1+0.1*delta).*cos(theta);
    subplot(1,2,2), plot([-x0,x0],[y0,y0],'k',[-x,x],[y,y],'r')
    axis('square')
    axis([-1.5 1.5 -1.5 1.5])
    title('Sphere Drainage Contour')
%    title(['Surface Contour at t = ',num2str(t)])
    drawnow
end
pause(5) %Wait 5 seconds and repeat movie
end
