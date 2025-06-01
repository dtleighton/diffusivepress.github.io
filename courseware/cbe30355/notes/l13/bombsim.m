%% Simulation of an explosion
% In this script we simulate the autocatalytic instability leading to a
% nuclear explosion.  Simplistically, we will take the neutrons to be point
% particles which are diffusing through a fissile material.  At each time
% step they take a random step (diffusion) in the x, y, and z directions.
% The length of this step is a normally distributed random variable with
% standard deviation (2*D*dt)^.5 where D is the diffusion coefficient.
% Because we've rendered everything dimensionless, the diffusion
% coefficient is just 1. They also have a probability of duplicating
% themselves if a random number is less than dt*alpha (which is a parameter
% much less than one). Particles which leave the sphere of dimensionless
% radius of one are eliminated from the list of particles.
%
% We shall use a value of alpha which is a bit greater than the critical
% value of pi^2 so that it blows up (eventually).  We start with our
% particles distributed uniformly over the sphere.  Because most of them
% are close to the edges, the total number initially decreases a bit (a
% uniform distribution is not the most unstable mode).  Eventually,
% however, we wind up with the distribution corresponding to the most
% unstable mode and the number of particles then increases exponentially in
% time.
%
% Note that while this is written for a specific application (a bomb),
% similar approaches can be used for many other autocatalytic problems:
% viral spreading is a good example!

% We kick off the simulation with a hundred neutrons.
n=100;
x=2*rand(n,1)-1;
y=2*rand(n,1)-1;
z=2*rand(n,1)-1;
r=(x.^2+y.^2+z.^2).^.5;

% We need to only keep the particles which are inside the sphere!
ikeep=find(r<1);
x=x(ikeep);
y=y(ikeep);
z=z(ikeep);

ntot=length(x); %The total number of particles inside the sphere

% We select a value of dt and alpha:

alpha = 15;
dt = 0.0002;

t = 0.0001; %We intialize time so that the graphics (title) looks better.

while t < 0.6

% We duplicate particles where a random number between 0 and 1 exceeds the
% threshold.  Note that this requires the threshold to be much less than 1.

dup = rand(size(x));
idup = find(dup<alpha*dt); %The ones which are duplicated
x=[x;x(idup)]; %We add these into our arrays
y=[y;y(idup)]; %note that matlab wants you to predimension the arrays, but
z=[z;z(idup)]; %this just makes it more complicated (but faster to run)...

% We add in the random walk:
x = x + (2*dt)^.5*randn(size(x)); %We do this for each direction
y = y + (2*dt)^.5*randn(size(x));
y = y + (2*dt)^.5*randn(size(x));

% We get rid of the particles which escape after the update:
r=(x.^2+y.^2+z.^2).^.5;
ikeep=find(r<1);
x=x(ikeep);
y=y(ikeep);
z=z(ikeep);

% We keep an array of the number of particles.  Again, matlab wants you to
% predimension the array (this is easier to do here, as you know how long
% it should be from the number of iterations), but we're not going to
% bother...
ntot=[ntot;length(x)];

% And we plot it up! (This is the time sink)

figure(1)
plot3(x,y,z,'o')
axis('square')
title(['Simulation for t = ',num2str(t),', N = ',num2str(length(x))])

% We want to add in a mesh sphere so we can see the domain more easily:
% (this bit of code was grabbed off the web by googling "3d plot of a 
% sphere in matlab".  That's an easy way to get tips on plotting things
% up!)
hold on
rp = ones(50, 50);
[th, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(-pi, pi, 50));
[xp,yp,zp] = sph2cart(th, phi, rp);
lightGrey = 0.8*[1 1 1]; % It looks better if the lines are lighter
surface(xp,yp,zp,'FaceColor', 'none','EdgeColor',lightGrey)
hold off

drawnow

t = t+dt; %We update time and do it again!

end

%% Growth rate
% We expect the number of particles to grow exponentially with a rate given
% by alpha-pi^2.  We can check this based on the final number of particles
% (extrapolating the exponential growth backwards in time).  It actually
% works pretty well.  Because it is a random simulation, however, it will
% look a bit different every time you run it!  It also behaves better for
% smaller values of dt (so that the dimensionless random walk length is
% smaller relative to the sphere size). You can run it much faster by
% commenting out the graphics in the first part, as it is the graphics
% which takes up all the computation time.  An alternative is to do the
% graphics only every ten or hundred iterations: it would speed things up
% by one to two orders of magnitude!

figure(2)
tall=[1:length(ntot)]'*dt;
semilogy(tall,ntot,tall,ntot(end)*exp((alpha-pi^2)*(tall-t)))
grid on
xlabel('time')
ylabel('number of particles')
legend('N','expected growth','Location','NorthWest')
title('Comparison of Growth Rate with Linear Stability Theory')

%% Distribution
% We can take the final distribution of particles and compare them to the
% expected most unstable mode.  The easiest way to do this is to sort the
% radial position of the particles and plot up the cumulative sum.  This is
% essentially the fraction of particles with radial locations less than r.
% For the most unstable mode calculated from linear stability theory this
% should match the function sin(pi*r)/pi-r*cos(pi*r).  As can be seen from
% the plot below this also matches quite well.

cumdist=[1:length(r)]'/length(r);
rpp=[0:.01:1];
expdist=sin(pi*rpp)/pi-rpp.*cos(pi*rpp);
figure(3)
plot(sort(r),cumdist,rpp,expdist)
axis([0 1 0 1])
xlabel('r')
ylabel('fraction particles less than r')
legend('MC simulation','Most Unstable Mode','Location','NorthWest')
title('Comparison of MC Simulation Distribution with Linear Stability Theory')

%% Conclusion:
% This example can be made to be much more realistic if you include effects
% such as partial reflection at the outer wall (used in making bombs,
% actually), and can be much faster to run if you predimension arrays and
% take advantage of symmetries in the system (you really only need to
% simulate one quadrant of the sphere) as well as eliminating the graphics.
% Monte Carlo simulations are very easy to write and with fast computers
% are a really convenient way to get quick answers to complicated problems.
% This is particularly true for non-linear problems for which there is no
% analytic solution, and for which other techniques such as finite
% difference or finite element calculations are more complex and demanding.