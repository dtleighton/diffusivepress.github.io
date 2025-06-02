%% Problem of the Day 08: Neumann Stability - Slab with heat generation
% In this problem we compare a finite difference marching solution to the
% exact result for a slab with uniform heat generation.  We are interested
% in the temperature at the bottom (insulated wall) as a function of time.
% The problem admits a nice closed form SL solution, so this gives us a
% point of comparison for the numerical result.  We look at time spacings
% just above and below the Neumann stability criterion.

%% The exact solution:
% We have the Sturm Liouville solution:

tsl = [0:.001:2];
n = [1:100]'; % we use 100 eigenvalues (overkill)
sigma = (n-.5)*pi; % the eigenvalues

Tbotsl = 1/2 + sum(2*(-1).^n./sigma.^3.*exp(-sigma.^2*tsl));

figure(1)
plot(tsl,Tbotsl)
xlabel('t')
ylabel('T|y = 0')
title('Temperature at the bottom')
grid on

%% The marching solution below the Neumann condition
% We use the center difference Euler method marching solution.  We choose a
% discretization of n in the spatial domain and thus have a time
% discretization of less than 0.5/n^2 for stability:

n = 20
dy = 1/n;

a = (diag(ones(n,1),-1)+diag(ones(n,1),1)-2*diag(ones(n+1,1)))/dy^2;

format short e
dt = 0.50*dy^2 %The maximum dt for stability

tkeep = [0:dt:3]; %The times we keep
Tbot = zeros(size(tkeep)); %We initialize the bottom temperature array

Tbot(1) = 0;
T = zeros(n+1,1); %our initial temperature distribution

for i = 2:length(tkeep)
    T = T + dt * (a * T + 1); %We add in the source
    T(n+1) = 0; %The upper BC
    T(1) = 4/3*T(2) - 1/3*T(3); %The lower BC
    Tbot(i) = T(1); %We keep the bottom temperature
end

figure(1)
hold on
plot(tkeep,Tbot,'--')
hold off
legend('Sturm-Liouville Solution','Marching Solution')
axis([0 2 0 .5])

%% Now we increase dt by just a bit...
dt = 0.51*dy^2 %just above stability!
format short

Tbotnew = zeros(size(tkeep)); %We initialize the bottom temperature array

Tbotnew(1) = 0;
T = zeros(n+1,1); %our initial temperature distribution

for i = 2:length(tkeep)
    T = T + dt * (a * T + 1); %We add in the source
    T(n+1) = 0; %The upper BC
    T(1) = 4/3*T(2) - 1/3*T(3); %The lower BC
    Tbotnew(i) = T(1); %We keep the bottom temperature
end

figure(1)
hold on
plot(tkeep,Tbotnew,':')
hold off
legend('Sturm-Liouville Solution','Marching Solution','Unstable Solution')
axis([0 2 0 .5])
