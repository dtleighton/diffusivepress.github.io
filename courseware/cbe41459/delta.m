function y=delta(x)
%This function integrates the function fderiv.m from zero to some value
%etalimit, and returns the deviation between the boundary conditions at
%etalimit and the correct values.  In order to run in the student version
%we return a single scalar deviation (e.g., fmins or fminsearch).  We can
%watch convergence to the solution graphically through the string variable
%"flag".  We use the canned integration routine ode23.m.

global etalimit flag

f0=[0 x(1) 0 x(2) 0]; %The initial condition with the two shooting parameters

[eta,f]=ode23('fderiv',[0 etalimit],f0);

if flag=='y'
 figure(1)
 plot(eta,f);%comment this out to run faster!
 legend('f1','f2','f3','f4','f5');
 title(['Plot for u(0) = ',num2str(x(1)),' T(0) = ',num2str(x(2))]);
drawnow;
end
nlast=length(eta);
y=f(nlast,2)^2+(f(nlast,5)-.5)^2; %The scalar "miss" that we are trying to minimize.

