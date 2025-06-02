%% The Nusselt-Graetz Problem: Constant Heat Flux at the Wall
% We solve the Sturm-Liouville problem for constant heat flux at the wall
% for laminar flow through a circular tube.
%
% The eigenvalue problem for the decaying solution is
%
% ((x*y')' + lambda*2*r*(1-r^2)*y = 0
%
% y'(0) = 0; y'(1) = 0
%
% The asymptotic solution at large z is:
%
% Tinf = 2*z + r.^2 - r.^4/4 - 7/24
%
% Thus, the initial value for the Sturm-Liouville expansion is just the
% negative of this function.
% 
% Solving the problem we get:

p = @(x) x;
q = @(x) zeros(size(x));
w = @(x) 2*x.*(1-x.^2);
bc = [0,1,0,1];

n = 100; %The number of points we would like (the number of intervals)

[lambda, eigenvecs] = slsolve(p,q,w,bc,n);

Tinf = @(z,r) 2*z + r.^2 - r.^4/4 - 7/24

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
    numerator = -weights*(w(r).*Tinf(0,r).*eigenvecs(:,i));
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

%% Temperature of the Cylinder Wall
% We want the temperature at the outer radius, necessary for obtaining the
% heat transfer coefficient:

z = [0.0005:.001:1];
twall = zeros(size(z)); %We initialize the array

for i = 1:length(z)
    twall(i) = Tinf(z(i),1) + sum(a.*exp(-lambda*z(i)).*eigenvecs(end,:)');
end

figure(2)
plot(z,twall)
xlabel('dimensionless axial position')
ylabel('cylinder wall temperature')
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
        tprofile(i,j) = Tinf(zplot(j),r(i)) + sum(a.*exp(-lambda*zplot(j)).*eigenvecs(i,:)');
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
% with axial position.  The dimensionless value is just the inverse of
% twall - Tb(z), where Tb(z) is just 2*z.  So:

h = 1./(twall - 2*z);

figure(4)
plot(z,h,z,24/11*ones(size(z)),'--')
axis([0,1,0,4])
xlabel('z')
ylabel('dimensionless heat transfer coefficient')
legend('transient solution','asymptotic solution')
grid on

%% Conclusion
% The initial heat transfer coefficient is larger than the asymptotic value
% due to the small thickness of the developing boundary layer. After a
% short distance, however, the asymptotic result is reached.  The zero
% initial eigenvalue is the result of a zero derivative condition at both r
% = 0 and r = 1, however the coefficient for it is likewise zero, removed
% via the energy balance enforced on the asymptotic solution.