%% Simulation of a Quenched Cube
% In this simulation we use Monte Carlo integration to determine the
% average temperature of a cube as a function of time.  Its initial
% temperature is 1 and the surface temperature is zero.  We start with N
% tracers distributed over the positive chunk of the cube, reflect at the
% symmetry planes, and remove tracers which exceed 1 in any direction.  The
% average temperature is just the number of tracers remaining in the
% domain!  Note that the step size in time must be small enough that the
% (2*dt)^1/2 jump is small with respect to the domain.

N = 10000; % We start with lots of tracers

dt = 0.00005; % This yields a random walk step of 1%

x = rand(N,1); % The initial values of x
y = rand(N,1);
z = rand(N,1);

tfinal = .5; % How long we run it for.

tall = [0:dt:tfinal];

t = 0;
i = 1; % A counter.

isout = max(max(x,y),z); %This returns the maximum of x, y, and z for each
ikeep = find(isout<1);
nleft = N;

Tavgkeep = zeros(size(tall));

Tavg = 1; % Our initial temperature

Tavgkeep(1) = Tavg;

while nleft>0 & t<tfinal
    i = i+1;
    t = tall(i);
    x(ikeep) = x(ikeep) + (2*dt)^.5*randn(nleft,1);
    y(ikeep) = y(ikeep) + (2*dt)^.5*randn(nleft,1);
    z(ikeep) = z(ikeep) + (2*dt)^.5*randn(nleft,1);
    x = abs(x); % We reflect at zero
    y = abs(y);
    z = abs(z);
    isout = max(max(x,y),z); %This returns the maximum of x, y, and z for each
    ikeep = find(isout<1);
    nleft = length(ikeep);
    Tavgkeep(i) = nleft/N;
end

figure(1)
plot(tall, Tavgkeep)
xlabel('t')
ylabel('Average Temperature')
grid on

