%% Illustration of Numerical Solution Techniques
% In this script we illustrate the numerical solution techniques of the
% Euler Method, the 2 stage Runge-Kutta method, the 4 stage R-K method, and
% an adaptive integration method built into matlab, ode23.  We choose as an
% example an undamped forced oscillator: a problem for which there is a
% simple analytic solution.  Our problem is:
%
% y'' + y = sin(a*t)
%
% y(0) = y'(0) = 0
%
% The behavior of this oscillator depends on a, the dimensionless driving
% frequency.  The amplitude blows up as a approaches 1, the natural
% frequency of the oscillator.
%
% We solve this problem as a pair of first order equations, such that y(1)
% is y and y(2) is y'.  The corresponding derivatives are:
%
% ydot = @(t,y) [y(2) ; -y(1) + sin(a*t)];
%
% With initial conditions: y0 = [0 ; 0];
% 
% We have the analytic solution to this equation:
%
% yexact = @(t) (sin(a*t) - a* sin(t))/(1-a^2);

%% The Euler Method
% We plot things up over four periods (e.g., up to 8pi).  We just define
% the derivative and the time step.  We also specify a, the frequency of
% the driving term.

a = 0.2;

% The exact solution:

yexact = @(t) (sin(a*t) - a* sin(t))/(1-a^2);

dt = 0.2;

tall = [0:dt:8*pi]';

y0 = zeros(2,1); % The initial condition

ydot = @(t,y) [y(2) ; -y(1) + sin(a*t)];

yem = zeros(2,length(tall)); % we keep both y1 and y2.

n = length(tall)-1; %the number of steps

for i = 1:n
    yem(:,i+1) = yem(:,i) + dt*ydot(tall(i),yem(:,i));
end

figure(1)
plot(tall,yexact(tall),'-o',tall,yem(1,:))
xlabel('t')
ylabel('y')
legend('exact solution','Euler Method')

%% The Two-Stage Runge-Kutta Technique
% Very little needs to be added to the Euler method for the 2s R-K
% technique.  The error, for this step size, is very small and is
% graphically indistinguishable from the exact solution.

y2s = zeros(2,length(tall)); % we keep both y1 and y2.
for i = 1:n
    k1 = dt*ydot(tall(i),y2s(:,i));
    k2 = dt*ydot(tall(i+1),y2s(:,i)+k1);
    y2s(:,i+1) = y2s(:,i) + (k1+k2)/2;
end

figure(1)
hold on
plot(tall,y2s(1,:),'r')
hold off
legend('exact solution','Euler Method','2 stage R-K')

%% The Four-Stage Runge-Kutta Technique
% This is very similar to the two stage technique, except there are four
% intermediate steps.  The method is even more accurate.

y4s = zeros(2,length(tall)); % we keep both y1 and y2.
for i = 1:n
    k1 = dt*ydot(tall(i),y4s(:,i));
    k2 = dt*ydot(tall(i)+dt/2,y4s(:,i)+k1/2);
    k3 = dt*ydot(tall(i)+dt/2,y4s(:,i)+k2/2);
    k4 = dt*ydot(tall(i)+dt,y4s(:,i)+k3);
    y4s(:,i+1) = y4s(:,i) + (k1+2*k2+2*k3+k4)/6;
end

figure(1)
hold on
plot(tall,y4s(1,:),'g')
hold off
legend('exact solution','Euler Method','2 stage R-K','4 stage R-K')

%% Adaptive Integrator
% We use the adaptive integration method ode23.m supplied with matlab.  It
% is also very accurate, and would require fewer steps.

[tout yout] = ode23(ydot,[0 8*pi],y0);

figure(1)
hold on
plot(tout,yout(:,1),'k')
hold off
legend('exact solution','Euler Method','2 stage R-K','4 stage R-K','ode23')

%% Error Comparison
% A better way of looking at the different methods is to determine the
% maximum error of each.  Here, since we have the exact solution, we can
% simply subtract it off and determine the maximum deviation.

errorem = max(abs(yem(1,:)-yexact(tall')))

error2s = max(abs(y2s(1,:)-yexact(tall')))

error4s = max(abs(y4s(1,:)-yexact(tall')))

errorode23 = max(abs(yout(:,1)-yexact(tout)))

%% Other Approaches
% Often it is desired to integrate to some condition (e.g., where a
% function value reaches zero) rather than over a discrete range in time.
% This can be easily done for RK techniques using a while loop rather than
% a for loop (updating the index at each step).  At the conclusion of the
% integration the array of function values and independent variables is
% trimmed (it needs to be predimensioned for efficiency) and the last
% element is adjusted via interpolation for improved accuracy.  Note that
% you can also adjust the array size inside the while loop as well.  This
% can also be done using the adaptive integrator ode23 via the "options"
% command.  An implementation of a "while" loop approach is given below.

t0 = 0; %The initial time
y0 = [0; 0]; %The inital value (column vector)

nchunk = 20; %We will predimension the arrays and expand as necessary
tallw = zeros(1, nchunk); %We will keep time as a row vector
y4sw = zeros(length(y0),nchunk); %The y values are an array

tallw(1) = t0;
y4sw(:,1) = y0;

i = 1;

while (y4sw(1,i)>0)||(i == 1); %The truncation condition, doing it at least once
    
    if i+1 > length(tallw); %Redimensioning the array
        y4sw = [y4sw, zeros(length(y0),nchunk)];
        tallw = [tallw,zeros(1, nchunk)];
    end
    
    k1 = dt*ydot(tallw(i),y4sw(:,i));
    k2 = dt*ydot(tallw(i)+dt/2,y4sw(:,i)+k1/2);
    k3 = dt*ydot(tallw(i)+dt/2,y4sw(:,i)+k2/2);
    k4 = dt*ydot(tallw(i)+dt,y4sw(:,i)+k3);
    
    y4sw(:,i+1) = y4sw(:,i) + (k1+2*k2+2*k3+k4)/6;
    tallw(i+1) = tallw(i) + dt;

    i = i+1; %We update i
end

tallw = tallw(1:i); %We trim the values
y4sw = y4sw(:,1:i);

% Now for the interpolation for the root crossing:
f = y4sw(1,end-1)/(y4sw(1,end-1)-y4sw(1,end)); %Test on y = 0
tallw(end) = tallw(end-1)+f*dt;
y4sw(:,end) = y4sw(:,end-1) + f*(y4sw(:,end)-y4sw(:,end-1));

figure(2)
plot(tall,yexact(tall),tallw,y4sw(1,:),'*')
xlabel('t')
ylabel('y')
legend('exact solution','4 stage R-K')
text(5,.5,['Root crossing at t = ',num2str(tallw(end))])
%% Conclusion
% As can be seen, the Euler Method isn't very accurate for this
% differential equation.  The two RK methods are much more accurate, with
% the 4 stage method being three orders of magnitude better than the 2
% stage method for this choice of dt.  The accuracy of the adaptive
% quadrature routine is determined by its tolerances, and can be adjusted.
% For the choice of parameters used here it has the same total number of
% steps as the other methods and lies between the 2 stage and 4 stage
% methods in accuracy.