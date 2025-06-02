%% Sturm-Liouville Example: Quenching of a Sphere at Finite Bi
% We examine the numerical solution to the Sturm-Liouville problem of a
% quenched sphere of radius 1 with finite Biot number.  The asymptotic
% temperature is zero (it has already been subtracted off), so the SL
% problem is:
%
% (x^2*y')' + lambda*x^2*y = 0
%
% y'(0) = 0; y(1) + 1/Bi * y'(1) = 0
%
% Solving the problem we get:

Bi = 2

p = @(x) x.^2;
q = @(x) zeros(size(x));
w = @(x) x.^2;
bc = [0,1,1,1/Bi];

n = 100; %The number of points we would like (the number of intervals)

[lambda, eigenvecs] = slsolve(p,q,w,bc,n);

% And that's it!

%% Eigenvalues, Eigenvectors, and Coefficients
% We are interested in the lead eigenvalues, coefficients, and
% eigenvectors.  We just look at the first five:

firsteigenvecs = lambda(1:5)

% And we calculate the coefficients using the Trapezoidal Rule:

x = [0:1/n:1]';
Tinit = ones(size(x));

weights = ones(1,n+1);
weights(1) = 0.5;
weights(n+1) = 0.5;
weights=weights/n;

a = zeros(length(lambda),1);

for i = 1:length(lambda)
    numerator = weights*(w(x).*Tinit.*eigenvecs(:,i));
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

%% Temperature of the Sphere Center
% We want the temperature at the center.  This will be Tinf + Tdecaying 
% evaluated at x = 0.  Tinf is just zero, so:
t = [0.0005:.001:2.5];
tcenter = zeros(size(t)); %We initialize the array

for i = 1:length(t)
    tcenter(i) = sum(a.*exp(-lambda*t(i)).*eigenvecs(1,:)');
end

figure(2)
plot(t,tcenter)
xlabel('dimensionless time')
ylabel('center temperature')
grid on

%% Temperature Profile at Various Times
% We can also plot up the temperature distribution for specific times.  You
% will note the issue near the origin at very short times.  This is known
% as the Gibbs ringing phenomenon and is well known in signal processing.

tplot = [0.0002,.05,.1,.2,.4,.8,1.6,10]';

tprofile = zeros(length(x),length(tplot));
for j = 1:length(tplot)
    for i=1:length(x)
        tprofile(i,j) = sum(a.*exp(-lambda*tplot(j)).*eigenvecs(i,:)');
    end
end
figure(3)
plot(x,tprofile)
legend(num2str(tplot),'Location','West')
xlabel('r')
ylabel('Temperature')
title('Temperature Distribution at Different Times')
grid on
        