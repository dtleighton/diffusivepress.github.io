%% Sturm-Liouville Example: Varying Heat Capacity and Conductivity
% We examine the numerical solution to the Sturm-Liouville problem of a
% slab of width 1 undergoing a uniform heat flux at the bottom, and a fixed
% temperature at the top.  The conductivity k and the heat capacity are
% linear functions of position.  The dimensionless Sturm-Liouville problem
% for the decaying solution is:
%
% ((1+c2*x)*y')' + lambda*(1+c1*x)*y = 0
%
% y'(0) = 0; y(1) = 0
%
% The asymptotic solution at long times is the logarithm:
%
% Tinf = (1/c2)*log((1+c2)./(1+c2*x))
%
% Thus, the initial value for the Sturm-Liouville expansion is just the
% negative of this function.
% 
% Solving the problem we get:

c1 = 2;
c2 = -.7;

p = @(x) (1 + c2*x);
q = @(x) zeros(size(x));
w = @(x) (1 + c1*x);
bc = [0,1,1,0];

n = 100; %The number of points we would like (the number of intervals)

[lambda, eigenvecs] = slsolve(p,q,w,bc,n);

Tinf = @(x) (1/c2)*log((1+c2)./(1+c2*x))

% And that's it!

%% Eigenvalues, Eigenvectors, and Coefficients
% We are interested in the lead eigenvalues, coefficients, and
% eigenvectors.  We just look at the first five:

firsteigenvecs = lambda(1:5)

% And we calculate the coefficients using the Trapezoidal Rule:

x = [0:1/n:1]';

% The Trapezoidal Rule weights:
weights = ones(1,n+1);
weights(1) = 0.5;
weights(n+1) = 0.5;
weights=weights/n;

a = zeros(length(lambda),1);

for i = 1:length(lambda)
    numerator = -weights*(w(x).*Tinf(x).*eigenvecs(:,i));
    denominator = weights*(w(x).*eigenvecs(:,i).^2);
    a(i) = numerator/denominator;
end

firstcoefficients = a(1:5)

% And we plot the first five eigenfunctions:
figure(1)
plot(x,eigenvecs(:,1:5))
xlabel('x')
ylabel('y')
title('First Five Eigenfunctions')
legend('n = 1','n = 2','n = 3','n = 4','n = 5')
grid on

%% Temperature of the Lower Wall
% We wanted the temperature of the lower wall.  This will be Tinf +
% Tdecaying evaluated at x = 0.  So:
t = [0.0005:.001:5];
tbottom = zeros(size(t)); %We initialize the array

for i = 1:length(t)
    tbottom(i) = log(1+c2)/c2 + sum(a.*exp(-lambda*t(i)).*eigenvecs(1,:)');
end

figure(2)
plot(t,tbottom)
xlabel('dimensionless time')
ylabel('bottom wall temperature')
grid on

%% Temperature Profile at Various Times
% We can also plot up the temperature distribution for specific times.  You
% will note the issue near the origin at very short times.  This is known
% as the Gibbs ringing phenomenon and is well known in signal processing.

tplot = [0.0002,.05,.1,.2,.4,.8,1.6,3.2]';

tprofile = zeros(length(x),length(tplot));
for j = 1:length(tplot)
    for i=1:length(x)
        tprofile(i,j) = log((1+c2)/(1+c2*x(i)))/c2 + sum(a.*exp(-lambda*tplot(j)).*eigenvecs(i,:)');
    end
end
figure(3)
plot(x,tprofile,x,Tinf(x))
legend(num2str(tplot),'asymptote')
xlabel('x')
ylabel('Temperature')
title('Temperature Distribution at Different Times')
grid on
        