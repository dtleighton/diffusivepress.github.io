%% Example of MC error propagation
% In this example we compare error propagation calculated via MC with that
% calculated via the standard (small variance) formula.  We apply it to the
% function f(x) = x1*log(x2).  The variables x1 and x2 are assumed to be
% independent normally distributed variables with known variance.

%% The standard formula
% This relies on a Taylor Series expansion of the function about expected
% values and requires the gradient of the function. We determine the 90%
% confidence interval (+/- 1.645 standard deviations).

x1 = 1;
x2 = 2;

sigx1 = 0.1;
sigx2 = 0.3;

f = @(x1,x2) x1.*log(x2); %The function
gradf = @(x1,x2) [log(x2),x1./x2]; %The gradient

varx = [sigx1^2,0;0,sigx2^2]; %The matrix of covariance of x

z = f(x1,x2)
varz = gradf(x1,x2)*varx*gradf(x1,x2)';
sigz = varz^.5
interval = z + sigz*[-1.645 1.645]

%% MC calculation
% We do the same calculation, but via MC simulation.  We add in a random
% perturbation to x1 and x2 and determine the resulting function values.
% We then sort this array and obtain the confidence interval.

n = 10000; %Lots of points...
x1mc = x1 + sigx1*randn(n,1);
x2mc = x2 + sigx2*randn(n,1);

zmc = f(x1mc,x2mc);
zmc = sort(zmc);

intervalmc = [zmc(n*.05),zmc(n*.95)]

%% Conclusion
% The results of the two calculations are close, with the difference (with
% such a large value of n) being because the standard deviations weren't
% very small, and thus the standard formula is a bit off.  They match
% perfectly if you use a smaller value of the standard deviation of the two
% variables.  Advantages of the MC approach are that you don't have to take
% the gradient, you can actually plot up the distribution, and that you can
% handle function values which are not normally distributed.  Disadvantages
% are the computation time, difficulty with dealing with covariance of the
% variables if they aren't independent (although this can be overcome),
% slightly different values every time you run it, and occasional errors if
% you feed the function with values which cause problems (e.g., values of
% x2 which are negative in this case due to a large standard deviation).