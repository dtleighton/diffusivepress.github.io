%% POD 22: MC Simulation of Taylor Dispersion
% We construct a simple Monte Carlo simulation of dispersion in a channel,
% first with no effect of the side walls, and then with the no-slip side
% walls causing a velocity reduction.  The effect of the side walls is
% approximated by causing all tracers within 5/8 of the edges to have a
% velocity of zero.  By symmetry, we shall only consider the top right
% quarter of the rectangular channel.  It has a dimensionless half-height
% of 1 and a dimensionless half width of ab (e.g., a/b).

n= 1000;
ab = 5;

% We distribute them randomly.
x = ab*rand(n,1);
y = rand(n,1);
z = zeros(n,1);

dt = 0.00005; %Our time step
tquit = 5;
t = [dt:dt:tquit]'; %Our times

uz = @(y) 1.5*(1-y.^2); %our velocity

varz = zeros(1/dt,1); %we keep track of the variance.

for i=1:tquit/dt-1
    dz1 = uz(y);
    
%    k = find(x>ab-5/8); dz1(k)=0; %The side wall fix...

    x = x + randn(n,1)*(2*dt)^.5;
    x = abs(x); %We reflect in the x-direction at the center
    x = min(x, 2*ab-x); %We reflect in from the right side
    
    y = y + randn(n,1)*(2*dt)^.5;
    y = abs(y);
    y = min(y, 2-y); %We reflect in from the top
    
    dz2 = uz(y);
    
%    k = find(x>ab-5/8); dz2(k)=0; %The side wall again.
    
    z = z + (dz1+dz2)/2*dt; %We update the velocity
    
    varz(i) = var(z);

end

    figure(1)
    subplot(2,2,1),plot(x,y,'o')
    xlabel('x')
    ylabel('y')
    title('x-y distribution')
    drawnow
    
    subplot(2,2,2),plot(x,z,'o')
    xlabel('x')
    ylabel('z')
    title('x-z distribution')
    drawnow
    
    subplot(2,2,3),plot(y,z,'o')
    xlabel('y')
    ylabel('z')
    title('y-z distribution')
    drawnow
    
    subplot(2,2,4),plot(t(1:i),varz(1:i),t(1:i),4/105*t(1:i))
    xlabel('t')
    ylabel('z variance')
    title('variance in the z direction')
    legend('MC simulation','TA dispersion calculation')
    drawnow

    
    
    