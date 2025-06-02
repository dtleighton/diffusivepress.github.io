%% Sturm-Liouville Example Problem
% We examine the numerical solution to the Sturm-Liouville problem of a
% slab of half-width 1 being quenched.  Initially the temperature is 1 and
% at t = 0 the surface at x = 1 is lowered to zero.  We use the symmetry
% condition at the centerline.  Thus:

p = @(x) ones(size(x));
q = @(x) zeros(size(x));
w = @(x) ones(size(x));
bc = [0,1,1,0];

n = 100; %The number of points we would like (the number of intervals)

[lambda, eigenvecs] = slsolve(p,q,w,bc,n);

% And that's it!

%% Comparison of Eigenvalues
% We have the exact eigenvalues ((n-1/2)*pi)^2.  We can compare them
% graphically and numerically.  We find that the first 11 eigenvalues match
% to within 1%, but then the deviation gets larger.  The highest
% eigenvalues are off by a factor of two or so.  This is generically the
% case: the eigenvalues (and eigenfunctions) are sensitive to the degree of
% discretization as each higher eigenfunction has one more zero crossing,
% and you lose numerical accuracy!  These higher eigenvalues would decay
% away very quickly, however.
figure(1)
neigenvec=length(lambda); 
exactlambda = (([1:neigenvec]'-.5)*pi).^2;
eratio = lambda./exactlambda;
plot([1:neigenvec],eratio,'o-')
xlabel('n')
ylabel('eigenvalue/exact eigenvalue')
grid on

first20ratios = eratio(1:20)

%% First Five Eigenfunctions
% We plot up the first five eigenfunctions.  You will note that the
% eigenfunctions have a number of zero crossings that increase by one each
% time.

x = [0:1/n:1]';
figure(2)
plot(x,eigenvecs(:,1:5))
xlabel('x')
ylabel('y')
title('First Five Eigenfunctions')
legend('n = 1','n = 2','n = 3','n = 4','n = 5')
grid on

%% Calculation of the Coefficients
% We make use of orthogonality to calculate the coefficients.  We can use
% the trapezoidal rule to do the integration.  Comparison with the exact
% solution shows that the first coefficients are very good, but after
% awhile the values are a little random.  These coefficients are very
% small, however.

weights = ones(1,n+1);
weights(1) = 0.5;
weights(n+1) = 0.5;
weights=weights/n;

tinit = ones(n+1,1); %The initial temperature

a = zeros(length(lambda),1);

for i = 1:length(lambda)
    numerator = weights*(w(x).*tinit.*eigenvecs(:,i));
    denominator = weights*(w(x).*eigenvecs(:,i).^2);
    a(i) = numerator/denominator;
end

aexact = -2*(-1).^[1:neigenvec]'./([1:neigenvec]'-.5)/pi;

figure(3)
plot([1:neigenvec]',a,'o',[1:neigenvec]',aexact,'*')
xlabel('n')
ylabel('coefficient')
grid on
    
%% Comparison of Decaying Solution
% Here we calculate the temperature at the centerline as a function of
% time.  We have both the exact solution as well as our numerical one.  We
% plot this up for a dimensionless time of 4.  The match is pretty much
% perfect except exactly at t = 0.  This error occurs because of the
% problems with all of the higher eigenvalues and eigenfunctions - which
% decay away almost instantly.  This related to the classic Gibbs ringing
% phenomenon, which is well known in signal processing.  We avoid this by
% starting our time series at a value slightly greater than zero.

t = [0.0005:.001:4];
cltemps = zeros(size(t));
cltempexact = zeros(size(t));

for i = 1:length(t)
    nrange = [1:1000]; %We use lots of eigenvalues for the exact expression
    cltempexact(i) = sum(-2*(-1).^nrange./(nrange-.5)/pi.*exp(-t(i)*(nrange-.5).^2*pi^2));
    
    %And now for the numerical value:
    cltemps(i) = sum(a.*exp(-lambda*t(i)).*eigenvecs(1,:)');
end

figure(4)
plot(t,cltemps,t,cltempexact,'--')
xlabel('dimensionless time')
ylabel('centerline temperature')
legend('numerical solution','exact solution')
grid on

figure(5)
plot(t,cltemps,t,cltempexact,'--')
xlabel('dimensionless time')
ylabel('centerline temperature')
legend('numerical solution','exact solution')
axis([0 1e-2 .99 1.01])
grid on