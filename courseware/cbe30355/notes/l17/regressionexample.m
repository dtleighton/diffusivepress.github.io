%% Demonstration of Linear Regression / Error Analysis
% In senior lab students study the velocity due to natural convection
% from a heated wire.  There are many aspects to this experiments, but one
% result is the magnitude of the centerline velocity as a function of
% height above the wire at constant power.  Theory predicts that this
% should obey a power law, going as height ^ (1/5).  The question you can
% use regression to answer is whether the data is consistent with this
% theory.
%% Data
% Let z be the height above the wire and u be the observed velocities, all
% in cgs units.  We have from their measurements:

z = [1.49 2.76 4.03 5.30]'; % Make it a column vector!
u = [0.486 0.554 0.545 0.547]';

% and we plot it up:

figure(1)
loglog(z,u,'*')
grid on
xlabel('z (cm)')
ylabel('u (cm/s)')

%% Regression
% In linear regression, the model must be linear in the modeling parameters
% (not necessarily a line!), so for a power law we transform the model by
% taking the natural log.  Our matrix A is the matrix of modeling function
% (evaluated at different z) and our "right hand side" vector b is the data
% observed at each z.  Our vector x is the fitting parameters associated
% with the modeling functions.  So:

A = [log(z),ones(size(z))];
b = log(u);

x = A\b % Solved in a least squares sense!

figure(1)
hold on
plot(z,exp(A*x),'k')
hold off
legend('data','power law fit','Location','SouthEast')

%% Error Analysis
% It is always more complicated to calculate the uncertainty in a quantity
% than to get it the first time!  Here we use the standard method for
% getting the uncertainty in the slope (e.g., x(1) from our model) and
% calculate the 95% confidence interval based on the assumptions that:
%
% 1) The model is perfect: it really is a power law
%
% 2) The error in the log of the velocities is random and independent
%
% 3) The data points all have the same variance
%
% Note that all of these assumptions are questionable, which is why
% statistics is so slippery!
%
% We do the error calculation by applying the normal equations (which would
% yield the same values of x!) but which are more convenient for the error
% calculation:

k = inv(A'*A)*A';
x = k*b % Same values as before!

[n m] = size(A)

% The residual:
r = A*x - b;

% The data variance:
varb = r'*r/(n-m);

% The matrix of covariance of x:
varx = k*varb*k'

% The standard deviation of the slope x(1):
sigslope = varx(1,1)^.5

%% The Confidence Interval
% For a large number of data points, the 95% confidence interval is just
% +/- 2 sigma.  Because the number of degrees of freedom is so small,
% however, the confidence interval is governed by the t-distribution rather
% than the normal distribution.  You can look this up in tables, but matlab
% has the built in function tinv (inverse of the t distribution).  So:

confint = x(1) + sigslope*tinv([0.025 0.975],n-m)

%% Conclusion
% If you had just taken the +/- 2 sigma bounds, the theoretical slope of
% 0.2 would lie outside the bounds by a little bit.  Because of the very
% small number of data points, however, and the significant scatter, the
% 95% confidence interval contains the theoretical prediction.  This is
% true even at the 87% confidence level.  Because of the broad tails of the
% t-distribution for small numbers of degrees of freedom, you really need
% to have a larger number of data points to get a high percentage
% confidence interval.

conflevel = 1-2*tcdf((x(1)-0.2)/sigslope,n-m)






