%% Problem of the Day 11: Radial Fluid Expansion
% In this script we plot up the data from the class demonstration and
% compare it to our model.  We expect that the radius of the fluid ring
% raised to the 8th power will grow linearly in time, with a
% proportionality constant that depends on the fluid viscosity, initial
% volume, weight of the glass plate, etc.

t = [9 35 94 238]; %Ring crossing times in seconds
r = [6 7 8 9]*0.85; %Ring locations in cm.

% Some parameters:

mu = 14; %viscosity of glycerin in poise at 20°C
m = 582; %weight of top plate in grams
v0 = 3; %fluid volume in ml

tc = pi^3*r(1)^8*mu/(v0^2*m*980) %the characteristic time

tstar = (t-t(1))/tc; %dimensionless times
rstar = r/r(1); %dimensionless radial positions

figure(1)
plot(tstar,rstar.^8,'o',tstar,1+8/3*tstar,'--',tstar,1+1.56*8/3*tstar,'k')
xlabel('tstar')
ylabel('rstar^8')
legend('data','theory','theory with viscosity at 25°C','Location','NorthWest')
grid on

%% Conclusion
% The data certainly matches a linear increase in R^8 with time, however
% the growth is a little higher than would be predicted.  It is highly
% likely that this discrepancy is due to the viscosity being just a bit
% less than expected.  The viscosity of glycerin is a very strong function
% of temperature.  The data matches that at 25°C exactly - a viscosity of 9
% poise rather than 14 poise at 20°C.  This is quite reasonable as the
% plate was sitting on a light box which was putting out a significant
% amount of heat...