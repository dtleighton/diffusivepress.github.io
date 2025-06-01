%% Ethylene glycol viscosity
% Here we analyze data for the viscosity of aqueous ethylene glycol solutions.
% We take the data via cut and paste from the engineering toolbox.com:
%
% https://www.engineeringtoolbox.com/ethylene-glycol-d_146.html
% 
% The first column is the temperature in F, the second in C, and the
% remaining entries are dynamic viscosities at different compositions.  The
% concentrations are in volume fraction of ethylene glycol.  The data is
% fit to a quadratic polynomial in 1/T (deg K).  Because the data is in
% tablular form, you can determine the fitting coefficients for all the
% concentrations via linear regression in a single step!
%
% Once we have the fitting coefficients we can use them to calculate things
% such as the temperature coefficient (1/µ dµ/dT) very easily.  Careful
% examination of the temperature coefficient graph suggests that there may
% be a problem with the data for 65% volume fraction at lower temperatures!

vf = [25	30	40	50	60	65	100];

data=[40	4.4	3	3.5	4.8	6.5	9	10.2	48
80	26.7	1.5	1.7	2.2	2.8	3.8	4.5	15.5
120	48.9	0.9	1	1.3	1.5	2	2.4	7
160	71.1	0.65	0.7	0.8	0.95	1.3	1.5	3.8
200	93.3	0.48	0.5	0.6	0.7	0.88	0.98	2.4];

degf = data(:,1);
degc = data(:,2);

visc = data(:,3:end);

% We work in degrees K:

degk = degc + 273.16;

% We fit to a polynomial in 1/T:

a = [ones(size(degk)),1.0./degk,1.0./degk.^2];

x = a\log(visc)

figure(1)
semilogy(degc,visc,'o',degc,exp(a*x))
xlabel('temperature (deg C)')
ylabel('viscosity (cp)')
legend('25%','30%','40%','50%','60%','65%','100%')
title('viscosity of ethylene glycol solutions for different volume fractions')
grid on
zoom on

tempcoef=[zeros(size(degk)),-1.0./degk.^2,-2.0./degk.^3]*x

figure(2)
plot(degc,tempcoef)
xlabel('temperature (deg C)')
ylabel('temperature coefficient (1/deg K)')
legend('25%','30%','40%','50%','60%','65%','100%','Location','SouthEast')
title('temperature coefficient of ethylene glycol solutions for different volume fractions')
grid on


