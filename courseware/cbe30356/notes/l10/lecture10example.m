%% Temperature Distribution in a Finite Slab with Periodic Heating
% In this script we plot up the temperature distribution in a finite
% thickness slab where we have oscillatory heating at the bottom and an
% insulation condition at the top.  We start by plotting the temperature
% distribution at a number of times:

i = (-1)^.5; %just in case i was used as an index elsewhere!
beta = 4;

% We require t to be a column vector and y to be a row vector.

T = @(t,y) imag(exp(i*t)*(cosh(i^.5*beta*y)-sinh(i^.5*beta)/cosh(i^.5*beta)*sinh(i^.5*beta*y)))

t = [0,pi/2,pi,3*pi/2,2*pi]';

y = [0:.01:1];

figure(1)
plot(y,T(t,y))
grid on
xlabel('y')
ylabel('T')
legend(num2str(t))
title(['Plot of Temperature Distribution for beta = ',num2str(beta),' at different times'])

%% Plot of Upper Surface Amplitude as a Function of Beta
% We can also plot up the amplitude of the temperature variation of the
% upper surface as a function of beta.

beta = [0:.01:10];

Tamp = abs(1.0./cosh(i^.5*beta));

figure(2)
plot(beta,Tamp)
xlabel('beta')
ylabel('Upper Surface Amplitude')
title('Plot of Upper Surface Amplitude for Varying Womersley Numbers')

%% Conclusion
% As can be seen, for small values of the Womersley number the oscillations
% at the upper surface are the same as the lower: the whole slab is at a
% uniform, time varying temperature.  For large Womersley number, however,
% the amplitude of the oscillations at the upper surface decreases
% exponentially.