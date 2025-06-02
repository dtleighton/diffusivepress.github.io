%% POD 23: MC Simulation of Taylor Dispersion in AF4
% We construct a simple Monte Carlo simulation of dispersion in a
% asymmetric flow field flow fractionation device.  The distribution is
% exponential in y, and there is a simple shear flow in the x direction.
% We add in a second species to see separation based on the diffusivity.

n= 1000;

Dastar = 1.5; % A second species with a different diffusivity.
% We distribute them randomly.

y = -log(rand(n,1));
x = zeros(n,1);

ya = -Dastar*log(rand(n,1));
xa = zeros(n,1);

dt = 0.001; %Our time step
tquit = 500;
t = [dt:dt:tquit]'; %Our times

ux = @(y) y; %our velocity profile

varx = zeros(1/dt,1); %we keep track of the variance.
varxa = zeros(1/dt,1); %we keep track of the variance.

for i=1:tquit/dt-1
    dx1 = ux(y);
    y = y - dt+ randn(n,1)*(2*dt)^.5; %Updating y
    y = abs(y); %Reflection from the accumulating wall
    dx2 = ux(y);
    x = x + (dx1+dx2)/2*dt; %We update the velocity
    
    varx(i) = var(x);
    
    % Now we repeat for the second type of particle.

    dxa1 = ux(ya);
    ya = ya - dt+ randn(n,1)*(2*Dastar*dt)^.5;
    ya = abs(ya);
    dxa2 = ux(ya);
    xa = xa + (dxa1+dxa2)/2*dt; %We update the velocity
    
    varxa(i) = var(xa);
    
    % We plot up the distributions every thousand iterations.
    if i/1000==round(i/1000)
      figure(1)
      plot(x,y,'or',xa,ya,'ob')
      xlabel('x')
      ylabel('y')
      title(['x-y distribution at t = ',num2str(t(i))])
      drawnow
    end
    
end

figure(2)
plot(t(1:i),varx(1:i),'r',t(1:i),varxa(1:i),'b')
xlabel('t')
ylabel('x variance')
title('variance in the x direction')
legend('unit diffusivity',['diffusivity = ',num2str(Dastar)])

figure(3)
histogram([x,xa])
xlabel('x')
ylabel('frequency')
title(['Particle Locations at t = ',num2str(tquit)])

%% Conclusion
% As can be seen, the variance of the distribution with the higher
% diffusivity is much larger (going as D^3).  This actually degrades the
% separation significantly requiring a greater simulation time for complete
% separation.  The scalings are correct, however.  There is also
% considerable skewness to the distributions resulting from the shear flow.

