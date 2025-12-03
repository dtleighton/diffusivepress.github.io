%This script calculates the velocity and temperature distributions due to
%a heated wire.  The script requires the two functions delta.m and fderiv.m.
%The integration procedure uses the shooting method to guess the two initial
%conditions: the velocity and temperature at the centerline.  Because of the
%stiffness of the equations, the sensitivity of the integration procedure to
%the initial guesses for the shooting parameters requires that we "sneak up"
%on the final values.  We do this by gradually increasing the upper limit
%of integration etalimit (passed in as a global parameter) until the solution
%converges.  Integration itself is done using the canned routine ode23.m.
%The shooting method optimizations are run using the canned routine
%fminsearch.m.

%The program produces both a graph and a data file.  The datafile is saved
%under the name ['pr',num2str(pr),'.txt'], thus there will be one generated
%for every value of prandtl number used. Note that to recover the actual
%dimensional velocities you must calculate the scaling factors Uc and Vc,
%and you must multiply the vertical velocity (u) by x*^0.2 and divide the 
%horizontal velocity (v) by x*^0.4.  The temperature is scaled by that at the
%centerline.

%In addition to the final results, it is also possible to display graphics
%associated with the convergence process.  Note that this greatly slows down
%the program, so only do this if you either have a fast computer or really want
%to see it!  This option is controlled by the global variable "flag".

global pr etalimit flag

pr=input('Enter the Prandtl number: ');

flag=input('Do you want convergence graphics (y/n)? ','s');

etalimit=1/pr^.5; %a good starting point

x0=[1 pr^.5]; %initial guesses for the centerline velocity and temperature

dx=x0; %What we use to determine convergence

while norm(dx)/norm(x0)>0.001
 etalimit=etalimit*1.5;
%x=fmins('delta',x0); %Use this one for older versions of matlab
 x=fminsearch('delta',x0); %Use this one for newer versions of matlab
 dx=x-x0;
 x0=x;
end

%And now we generate the x and y velocities:
f0=[0 x(1) 0 x(2) 0]';
[eta f]=ode23('fderiv',[0 8],f0);

u=f(:,2); %The vertical velocity
v=-0.6*f(:,1)+0.4*eta.*f(:,2); %The horizontal velocity (antisymmetric in eta!)
temp=f(:,4); %The temperature

figure(1)
plot(eta,u,eta,v,'--',eta,temp/temp(1),'-.')
xlabel('eta')
ylabel('u* / x*^0^.^2, v* ^. x*^0^.^4, T/T(0)')
legend('u','v','T/T(0)')
title(['Plot of Dimensionless Velocity and Temperature for a Heated Wire, Pr = ',num2str(pr)])
text(3.1,.52,['dimensionless scaled centerline velocity = ',num2str(u(1))]);
text(3.1,.45,['dimensionless scaled centerline temperature = ',num2str(temp(1))]);
grid

data=[eta,u,v,temp];
savename=['pr',num2str(round(pr)),'.txt'];
save(savename,'data','-ascii')

% We use the self-similar solution to determine the ratio of stratification
% induced buoyancy to that from the heated wire.  We will do the integral
% of the dimensionless temperature exactly, rather than using the large Pr
% approximation discussed in the writeup.  The integrals are done using the
% mid-point rule, easily accurate enough for our purposes here.  Note that
% the tank temperature gradient at which stratification effects would be
% observed is Tc/Lx x*^-8/5 divided by the value computed below, which is
% about 3.

ig=(1-(f(:,1)/f(length(eta),1)).^(5/3));
deta=eta(2:length(eta))-eta(1:length(eta)-1);
igt=(ig(2:length(eta))+ig(1:length(eta)-1))/2;
igb=(temp(2:length(eta))+temp(1:length(eta)-1))/2;
stratratio=sum(igt.*deta)/sum(igb.*deta) %This is the ratio determining the allowable vertical stratification.
