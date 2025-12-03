% This script calculates the degree of subcooling of water necessary for
% the glycerin-water solution to come out at room temperature after mixing.
% It is based on the volume fraction of glycerin in the solution (rather
% than mass fraction) so the function density.m is required to convert
% things. It uses the function temprise.m to get the heat of mixing, and
% the function htcapacity.m to get the heat capacity of the solution.  It
% assumes the default temperature of 22 deg C.  The output is graphical for
% an array of volume fractions, but you can easily modify it to get the
% degree of subcooling for any desired volume fraction.  You can also adapt
% it to determine the amount of tap water and ice water to achieve the
% necessary degree of subcooling.

% OK, we start with volume fraction:

vf = [0:.05:.67]';

bf = 0.997; %a guess as to the true wt% of the glycerin (USP grade)

conc = vf*bf*density(bf)./(vf*density(bf)+(1-vf)*density(0));

deltacool = temprise(conc).*htcapacity(conc)/htcapacity(0)./(1-conc);
figure(1)
plot(vf,deltacool)
xlabel('volume fraction glycerin')
ylabel('required water subcooling')
title('Decrease of water temperature from ambient to balance heat of mixing')
grid on
zoom on

