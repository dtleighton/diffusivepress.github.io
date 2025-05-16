%This script produces a movie of startup flow in a pipe.  It uses the Euler
%method marching forward in time.

n=20; %The spatial discretization.
dr=1/n;
r=[0:dr:1];
u=zeros(size(r));
udot=zeros(size(r));

i=[2:n]; %The interior nodes
dt=0.5*dr^2; %We choose a time discretization for stability.

t=0;

while t<.5
 t=t+dt;
 
 % The derivative for the interior nodes
 udot(i)=1+((u(i+1)-u(i))/dr.*(r(i)+dr/2)...
    -(u(i)-u(i-1))/dr.*(r(i)-dr/2))/dr./r(i);

 u(i)=u(i)+udot(i)*dt; %We update interior nodes
 
 u(n+1)=0; % The velocity at the outer wall remains zero.

 % We use a zero derivative condition at the center.
 u(1)=2/3*(2*u(2)-0.5*u(3));

 figure(1)
 plot(u,r,'x',0.25*(1-r.^2),r,[0 0],[0 1])
 axis([0 .3 0 1])
 xlabel('u*')
 ylabel('r*')
 title(['Velocity Profile in a Pipe for t* = ',num2str(t)])
%  drawnow
 
 %Let's compare this to the exact solution!  We have the eigenvalues (roots
 %of J0, we keep ten):
sigma =[2.40482555769577
   5.52007811028631
   8.65372791291101
  11.79153443901428
  14.93091770848779
  18.07106396791092
  21.21163662987925
  24.35247153074930
  27.49347913204025
  30.63460646843198]';
 %Which yields the coefficients:
 an=-2.0./sigma.^3.0./besselj(1,sigma);
 %and the solution:
 uexact=0.25*(1-r.^2)+(an.*exp(-sigma.^2*t))*besselj(0,sigma'*r);
 %and we plot it up:
 hold on
 plot(uexact,r,'g')
 hold off
 
 %We can also plot up a movie of the deviation:
 figure(2)
 plot(r,u-uexact)
 xlabel('r*')
 ylabel('deviation')
 title(['Deviation from exact solution (ten eigenvalues) for t* = ',num2str(t)])
 drawnow
end

%Note that most of this stuff was the comparison...