%% Determination of Melt Location as a Function of Time
% Here we extend the problem of the day to determine the location at which
% our semi-infinite slab reaches a particular temperature.  Since we have
% the self-similar solution in terms of eta, it easiest to pick a
% particular time and plot the dimensional temperature as a function of
% y.  This is done below.  Note that we can pick tc to be anything we want.
% It is convenient to measure this in minutes.  We calculate the
% characteristic length using the thermal diffusivity of aluminum, and get
% the characteristic scaling for the temperature from the rated heat output
% of the source and the properties of the rod.

x = 1; % Our initial guess for the dimensionless temperature at the wall

x = fzero('miss2',x); % Our solution!

fdot = @(eta,f) [f(2); 0.5*(f(1)-eta*f(2))]; %The differential equation

f0 = [x,-1]; %The initial value

[etaout fout] = ode45(fdot,[0 10],f0);

% Now for the range of times.

t = [0.5:.5:4]';

y = t.^.5*etaout'; % This produces a matrix of values of y

T = t.^.5*fout(:,1)'; % A matrix of temperatures

% Now we introduce some properties
k = 152; %W/mûK thermal conductivity of aluminum 6061

alpha = 0.64; %cm^2/s thermal diffusivity

a = 0.00635; %m rod radius

Q = 40; %W heat source

delta = (alpha*60)^.5; %characteristic length scale in cm

q0 = Q/pi/a^2; %heat flux in W/m^2

Tc = q0*delta/100/k; %characteristic temperature rise (ûK)

T0 = 22; %ref temperature in ûC

T = T0 + T*Tc; %We convert the temperature to ûC

y = y*delta; %We convert the length to cm.

figure(1)
plot(y',T',[0,20],49*[1,1],'k--',[0,20],64*[1,1],'k--',[0,20],56.5*[1,1],'r--')
xlabel('y (cm)')
ylabel('T (ûC)')
legend(str2mat(num2str(t),'melt boundary','melt boundary','median'),'Location','bestoutside')
title('Dimensional Temperature Profile at Varying Times (min)')
grid on
axis([0 20 20 80])

%% Comparison to Sturm-Liouville Solution
% The self-similar solution is accurate for an infinite rod - but of course
% our rod isn't of infinite extent.  Thus, it will break down as the
% temperature propagates to the region near the end of the rod.
% Physically, this corresponds to the melting occurring much faster near
% the end of the experiment, as there is no longer a lot of "cold rod"
% beyond the end to absorb energy.
%
% There are two ways to deal with this.  The first is to simply add in a
% self-similar solution originating from a distance 2L from the original
% origin.  This woulld enforce the zero derivative (no flux) condition at
% the actual end of the rod, and would be good until the wave propagates to
% the far end (e.g., 2L).  A better way, however, is to solve the problem
% via asymptotic solution/separation of variables/Sturm-Liouville
% expansion.  The temperature profiles from these two methods at a few
% times are compared below.  As can be seen, agreement is very good for
% short times (as expected), but deviates significantly as the heat reaches
% the reflective boundary.

L = 19; %cm (rod length)
ystar = (L - y)/L; %dimensionless length based on the rod end (symmetry)
tc = L^2/alpha/60; %characteristic time in minutes
tstar = t/tc; %dimensionless time based on rod length
TcSL = Tc*L/delta; %conversion of Tc to our usual scaling

%We are going to need to use t as a scalar to get the function for Tstar
n = [1:100]'; %We use 100 eigenvalues (overkill)
Tstar = @(tt,y) tt + .5*y.^2 - 1/6 + sum((-2*(-1).^n/pi^2./n.^2.0.*exp(-n.^2*pi^2*tt)*ones(size(y))).*cos(n*pi*y))

figure(2)
nplotall = [2,4,6,8]';
for i=1:length(nplotall)
    nplot = nplotall(i);
    plot(y(nplot,:),T(nplot,:))
    hold on
end
legend(num2str(nplotall))

for i=1:length(nplotall)
    nplot = nplotall(i);
    plot(y(nplot,:),Tstar(tstar(nplot),ystar(nplot,:))*TcSL+T0,'k--')
    hold on
end
plot([0,L],56.5*[1,1],'r')
hold off
grid on
xlabel('y (cm)')
ylabel('T (ûC)')
legend(str2mat(num2str([t(nplotall);t(nplotall)]),'median melting point'),'Location','bestoutside')
title('Dimensional Temperature Profile at Varying Times (min)')
axis([0 L 20 80])

%% Conclusion
% Comparison with our crude measurements show that the melt times are
% longer than what would be obtained from the transient conduction model.
% The discrepancy may well be due to heat losses: the rod isn't
% particularaly well insulated!  Also, a fair amount of the power
% dissipated by the source isn't actually going into the rod, but instead
% is lost short of the connection.  Both effects would increase the time
% for melting to occur.
