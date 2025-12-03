% This script calculates the the sedimentation velocity of the tracer
% particles as a function of glycerin volume fraction

% OK, we start with volume fraction:

vf = [0:.05:.67]';

bf = 0.997; %a guess as to the true wt% of the glycerin (USP grade)

conc = vf*bf*density(bf)./(vf*density(bf)+(1-vf)*density(0));

mu = viscosity(conc);
rho = density(conc);
a = 0.0010; %particle radius

velocity = 2/9*(2.71-rho)*980*a^2.0./mu;

figure(1)
plot(vf,velocity)
xlabel('volume fraction glycerin')
ylabel('sedimentation velocity of tracers (cm/s)')
title('Sedimentation velocity of chalk particles')
grid on
zoom on

