%% The Nusselt-Graetz Problem: Constant Temperature at the Wall
% We solve the Sturm-Liouville problem for constant temperature at the wall
% for laminar flow through a circular tube.
%
% The eigenvalue problem for the decaying solution is
%
% ((x*y')' + lambda*2*r*(1-r^2)*y = 0
%
% y'(0) = 0; y(1) = 0
%
% The asymptotic solution at large z is just Tinf = 1.  Thus, the initial
% value for the Sturm-Liouville expansion is just -1.
% 
% Solving the problem we get:

p = @(x) x;
q = @(x) zeros(size(x));
w = @(x) 2*x.*(1-x.^2);
bc = [0,1,1,0];

n = 100; %The number of points we would like (the number of intervals)

[lambda, eigenvecs] = slsolve(p,q,w,bc,n);

Tinf = 1

% And that's it!

%% Eigenvalues, Eigenvectors, and Coefficients
% We are interested in the lead eigenvalues, coefficients, and
% eigenvectors.  We just look at the first five:

firsteigenvecs = lambda(1:5)

% And we calculate the coefficients using the Trapezoidal Rule:

r = [0:1/n:1]';

% The Trapezoidal Rule weights:
weights = ones(1,n+1);
weights(1) = 0.5;
weights(n+1) = 0.5;
weights=weights/n;

a = zeros(length(lambda),1);

for i = 1:length(lambda)
    numerator = -weights*(w(r).*Tinf.*eigenvecs(:,i));
    denominator = weights*(w(r).*eigenvecs(:,i).^2);
    a(i) = numerator/denominator;
end

firstcoefficients = a(1:5)

% And we plot the first five eigenfunctions:
figure(1)
plot(r,eigenvecs(:,1:5))
xlabel('x')
ylabel('y')
title('First Five Eigenfunctions')
legend('n = 1','n = 2','n = 3','n = 4','n = 5')
grid on

%% Flow Averaged Temperature and Heat Flux
% To obtain the heat transfer coefficient we need to know the flow averaged
% temperature and the heat flux at the wall.

z = [0.0005:.001:2];
Tavg = zeros(size(z)); %We initialize the array
q = zeros(size(z));
T = zeros(size(r));
dr = 1/n;
for j = 1:length(z)
    for i=1:length(r)
        T(i) = Tinf + sum(a.*exp(-lambda*z(j)).*eigenvecs(i,:)');
    end
    Tavg(j) = weights*(2*(1-r.^2).*2.*r.*T);
    q(j) = (1.5*T(end)-2*T(end-1)+0.5*T(end-2))/dr;
end

figure(2)
plot(z,Tavg,z,q)
xlabel('dimensionless axial position')
ylabel('flow averaged temperature, wall derivative')
legend('flow averaged temperature','wall derivative')
axis([0 2 0 2])
grid on

%% Temperature Profile at Various Times
% We can also plot up the temperature distribution for specific z
% locations.  You will note the issue near the origin at very small z.
% This is known as the Gibbs ringing phenomenon and is well known in signal
% processing.

zplot = [0.0002,.05,.1,.2,.4,.8]';

tprofile = zeros(length(r),length(zplot));
for j = 1:length(zplot)
    for i=1:length(r)
        tprofile(i,j) = Tinf + sum(a.*exp(-lambda*zplot(j)).*eigenvecs(i,:)');
    end
end
figure(3)
plot(r,tprofile)
legend(num2str(zplot))
xlabel('z')
ylabel('Temperature')
title('Temperature Distribution at Different Axial Positions')
grid on

%% Heat Transfer Coefficient
% Finally, we are interested in how the heat transfer coefficient varies
% with axial position.  The dimensionless value is the ratio of the heat
% flux to the difference between the bulk and wall temperatures.  So:

h = q./(1 - Tavg);
Nu = 2*h; %Convert to conventional definition hD/k

zstar = z/4; %Convert to conventional z/(UD^2/alpha)

figure(4)
plot(zstar,Nu,zstar,48/11*ones(size(zstar)),'--k')
axis([0,.1,2,7])
xlabel('z/(UD^2/alpha)')
ylabel('Nu = hD/k')
legend('transient solution','constant flux asymptotic solution')
grid on

%% Conclusion
% The initial heat transfer coefficient is larger than the asymptotic value
% due to the small thickness of the developing boundary layer. After a
% short distance, however, the asymptotic result is reached.  Both the heat
% transfer coefficient and the difference between wall and bulk temperature
% vanish, but the ratio is constant (and could be determined just from the
% lead eigenfunction).  The heat transfer coefficient is slightly less
% than that obtained from the constant wall heat flux boundary condition.
% This is the usual value used in correlations for the Nusselt number.