%% Separation of Variables Example: Startup of flow in a tube
% We examine the numerical solution to the startup of a fluid flowing
% through a circular tube: what happens in a vertical straw filled with
% fluid when you take your finger off the end.

% We have the asymptotic solution uinf:

uinf = @(r) (1-r.^2)/4;

% And we have the functions p, q, and w from the resulting spatial
% Sturm-Liouville problem:

p = @(x) x;
q = @(x) zeros(size(x));
w = @(x) x;

% The boundary conditions are zero derivative at the center and zero
% magnitude at r = 1:
bc = [0,1,1,0];

% We set the number of eigenvalues we would like (only the earlier ones are
% accurate) and the degree of discretization.

n = 100; %The number of points we would like (the number of intervals)

% and we use our solver:

[lambda, eigenvecs] = slsolve(p,q,w,bc,n);

% And that's it!

%% Eigenvalues, Eigenvectors, and Coefficients
% We are interested in the lead eigenvalues, coefficients, and
% eigenvectors.  We just look at the first five:

firsteigenvals = lambda(1:5)

% And we calculate the coefficients using the Trapezoidal Rule:

r = [0:1/n:1]';

% The Trapezoidal Rule weights:
weights = ones(1,n+1);
weights(1) = 0.5;
weights(n+1) = 0.5;
weights=weights/n;

a = zeros(length(lambda),1);

for i = 1:length(lambda)
    numerator = -weights*(w(r).*uinf(r).*eigenvecs(:,i));
    denominator = weights*(w(r).*eigenvecs(:,i).^2);
    a(i) = numerator/denominator;
end

firstcoefficients = a(1:5)

% And we plot the first five eigenfunctions:
figure(1)
plot(r,eigenvecs(:,1:5))
xlabel('r')
ylabel('y')
title('First Five Eigenfunctions')
legend('n = 1','n = 2','n = 3','n = 4','n = 5')
grid on

%% Velocity at the Centerline
% We are most interested in the velocity at the centerline.  This will be
% uinf + udecaying evaluated at r = 0.  So:
t = [0.0005:.001:1];
ucenter = zeros(size(t)); %We initialize the array

for i = 1:length(t)
    ucenter(i) = uinf(0) + sum(a.*exp(-lambda*t(i)).*eigenvecs(1,:)');
end

% We can also look at the one-term approximation:
ucenter1term = uinf(0) + a(1)*exp(-lambda(1)*t)*eigenvecs(1,1);

figure(2)
plot(t,ucenter,[0,.25],[0,.25],t,ucenter1term)
xlabel('dimensionless time')
ylabel('centerline velocity')
legend('numerical solution','short time asymptote','one term approx')
grid on

%% Velocity Profile at Various Times
% We can also plot up the velocity distribution for specific times.  You
% will note the issue near the origin at very short times.  This is known
% as the Gibbs ringing phenomenon and is well known in signal processing.

tplot = [0.0002,.05,.1,.2,.4,.8,1.6,3.2]';

uprofile = zeros(length(r),length(tplot));
for j = 1:length(tplot)
    for i=1:length(r)
        uprofile(i,j) = uinf(r(i)) + sum(a.*exp(-lambda*tplot(j)).*eigenvecs(i,:)');
    end
end
figure(3)
plot(r,uprofile,r,uinf(r))
legend(num2str(tplot),'asymptote')
xlabel('r')
ylabel('Velocity')
title('Velocity Distribution at Different Times')
grid on
        