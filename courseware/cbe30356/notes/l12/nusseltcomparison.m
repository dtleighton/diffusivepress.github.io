%% Comparison of Local Nusselt Number Calculations
% In this script (which must be run after the nusseltgraetz.m script)
% compares the asymptotic 48/11 solution, the Sturm-Liouville solution, and
% the thermal boundary layer solution for the local Nusselt number for the
% constant wall heat flux problem.  We convert z and h to the conventional
% scalings Nu = hD/k and z* = z/(UD^2/alpha).

zstar = z/4;

Nu = 2*h;

NuBL = 2*gamma(2/3)/9^(1/3)./zstar.^(1/3);

figure(1)
plot(zstar,Nu,zstar,48/11*ones(size(zstar)),zstar,NuBL,'--k')
xlabel('z/(UD^2/alpha)')
ylabel('Nu (local)')
legend('Exact Solution','Asymptotic Solution','Boundary Layer Solution')
title('Comparison of Constant Wall Heat Flux Nusselt Numbers')
axis([0 .08 3 10])
grid on

%% Conclusion
% Examination of the graph indicates that the boundary layer solution is
% good up to a dimensionless value of about 0.01 while the asymptotic
% solution is approached at a dimensionless value of 0.05.