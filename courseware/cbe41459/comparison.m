clear
format compact
% This script imports and plots up the velocity profile map for the
% sample pictures from the heated wire experiment.  It also compares 
% it to the theoretical profile.  Note that it is necessary for the
% data file produced from the theoretical program to be in this directory.

%We have the base name of the data files:
basename='sample'

%We load in the velocity data set:
name=[basename,'vel.mat'];
load(name,'-mat')

% We plot up the identified particle locations with symbols.  For the
% quiver plot we shift things back to the original (laboratory) reference
% frame.  We also move the lateral center to the peak velocity at the
% middle height.
ycshift = yc+[zeros(size(xc)),xc-median(xc)]*xfit;
drift = xfit(2) %The slope of the drift of the centerline

figure(1)
plot(ycshift,xc,'ok')
axis('equal')
ylabel('x (cm)')
xlabel('y (cm)')
title('Velocity Profile Quiver Plot')
labels=name; %We set up a list of labels.

%Now we add in the quiver velocity plot:
figure(1)
hold on
quiver(ycshift,xc,v,u,'.')
hold off

% We can also plot up the dimensional velocities as a function of y.
% Because it's also a function of x, we restrict ourselves to values of x
% close to the median of the data set.  Note that there is -a lot- more
% data in the file, so we can use the quiver plot (figure 1) to see what
% other values of x we can look at!  You will want to play with both xmed
% (the center of the data you are looking at for figures 2 and 3) and xd
% (the half-width of the vertical region you are looking at) to see how the
% profile varies with height for this data set.  You may even want to have
% xmed be negative (e.g., below the wire), although the theory isn't valid
% in that region, and only the horizontal velocity is of interest in that
% region.  Note that figures 4 and 5 are the dimensionless, scaled data and
% theory: they include all the data more than 1mm above the wire.

xmed=round(median(xc*10))/10 %We round off to the nearest mm.
xd=.2;
ii=find(abs(xc-xmed)<xd);
figure(2)
plot(yc(ii),u(ii),'o')
xlabel('y, cm')
ylabel('u, cm/s')
title(['Vertical Velocity Profile at x = ',num2str(xmed),'+- ',num2str(xd),'cm'])

figure(3)
plot(yc(ii),v(ii),'o')
xlabel('y, cm')
ylabel('v, cm/s')
title(['Horizontal Velocity Profile at x = ',num2str(xmed),'+- ',num2str(xd),'cm'])
grid on

% We also want to compare the dimensionless velocity profiles to theory.
% We have the characteristic velocities and length scales (taking Lx to
% be 1cm) for a particular experiment - you MUST calculate them for yours!:
Uc=.5142;
Ly=.2436;
Vc=Uc*Ly/1.0;
Tc=.6154;

ys=yc./xc.^0.4/Ly;
us=u./xc.^0.2/Uc;
vs=v.*xc.^0.4/Vc;

% We now plot up the scaled velocity profile:
kk=find(xc>.1); %We restrict these to partilces more than 1mm above the wire
figure(4)
plot(ys(kk),us(kk),'o')
xlabel('scaled dimensionless lateral position')
ylabel('scaled dimensionless vertical velocity')
title('Dimensionless Scaled Vertical Velocity')
axis([-10 10 -.25 1.25])
grid on

figure(5)
plot(ys(kk),vs(kk),'o')
xlabel('scaled dimensionless lateral position')
ylabel('scaled dimensionless horizontal velocity')
title('Dimensionless Scaled Horizontal Velocity')
labels2='Data';

% We can also compare to the theoretical velocity profile.  The data for
% This is stored in the file prXX.txt where XX is the prantl number:

pr=26;
name=['pr',num2str(pr)];
load([name,'.txt'])
data=eval(name);
eta=data(:,1);
ut=data(:,2);
vt=data(:,3);
eta=[-eta;eta];
ut=[ut;ut];
vt=[-vt;vt];
[eta,i]=sort(eta);
ut=ut(i);
vt=vt(i);
labels2=str2mat(labels2,'Theory');

%We add this profile to the dimensionless plots:
figure(4)
hold on
% We can introduce "scale factors" which scale the theoretical velocity and
% width by an arbitrary amount.  This is useful in determining the
% centerline velocity of the data which is often difficult to get if we
% don't have a measured velocity "right there".
yscale = 1; %These should be changed to match your data
uscale = 1;
plot(eta,ut,'k',eta*yscale,ut*uscale,'r')
legend(char(labels2,'Scaled Theory'))
text(2,.6,char(['velocity scale = ',num2str(uscale)],['width scale = ',num2str(yscale)]))
hold off

figure(5)
hold on
plot(eta,vt,'k')
legend(labels2)
hold off

%And we can add the profile to the dimensional plots as well:
etadim=eta*Ly*xmed^0.4;
udim=ut*Uc*xmed^.2;
vdim=vt*Vc/xmed^.4;

figure(2)
hold on
plot(etadim,udim)
plot(etadim*((xmed+xd)/xmed)^.4,udim*((xmed+xd)/xmed)^.2,':')
plot(etadim*((xmed-xd)/xmed)^.4,udim*((xmed-xd)/xmed)^.2,':')
plot(etadim*yscale,udim*uscale,'b--')
hold off

figure(3)
hold on
plot(etadim,vdim)
plot(etadim*((xmed+xd)/xmed)^.4,vdim/((xmed+xd)/xmed)^.4,':')
plot(etadim*((xmed-xd)/xmed)^.4,vdim/((xmed-xd)/xmed)^.4,':')
hold off

%We also estimate the slope of the horizontal velocity near the origin:
ikeep=find((abs(yc)<2*Ly*abs(xmed)^.4)&(xc>xmed-xd)&(xc<xmed+xd));
amat=[yc(ikeep),ones(size(yc(ikeep)))];
vfit=amat\v(ikeep);
dvdy=vfit(1)
dvdytheory=-Vc/Ly/xmed^.8*.23
figure(3)
hold on
plot(3*Ly*[-1;1],[3*Ly*[-1;1],[1;1]]*vfit,'k')
hold off
legend(char(labels2,['x = ',num2str(xmed+xd)],['x = ',num2str(xmed-xd)],'linear fit'))

% Now we turn to the tricky part: estimating the centerline velocity and
% profile width at a particular height above the wire.  This is often hard
% to do because you may not have a velocity measurement "right there", and
% you may have an incomplete profile at a particular height.  The procedure
% outlined below works pretty well -sometimes- and you always need to do a
% reality check to make sure it is giving you a reasonable value!

% Because the profile grows in width and magnitude as we move upwards, and
% the width of the profile "strip" we are looking at is finite, we will
% use the self-similar profile to adjust the vertical velocities and
% lateral positions to the target height:

uadj = u./(xc/xmed).^.2;
yadj = yc./(xc/xmed).^.4;

% Getting the centerline velocity at a particular height is actually pretty
% tricky.  This is because you have no direct control over where you are
% making the measurements (this occurs only where there happen to be
% particles) and because the profile is quite narrow.  We shall get an
% estimate by fitting a parabola and a triangle to the data near the
% centerline (we go out to a value of eta of about 0.6).  The parabola will
% underestimate the centerline velocity, and the triangle will tend to
% overestimate the velocity.  Comparing to the theoretical profile, it
% seems that a weighted average of these two values will give the best
% estimate.  Note that any spurious points will tend to distort the
% estimate, and that you have problems if the center of the profile isn't
% exactly at zero or if you have an insufficient number of points.  It is
% important to do a "reality check"!

ylim = 0.6*Ly*xmed^.4;
ikeep = find(abs(xc-xmed)<xd&abs(yadj)<ylim);
amat1 = [ones(size(yadj(ikeep))),yadj(ikeep).^2];
amat2 = [ones(size(yadj(ikeep))),abs(yadj(ikeep))];
uf1 = amat1\uadj(ikeep);
uf2 = amat2\uadj(ikeep);
ucl = 0.75*uf1(1)+0.25*uf2(1)
% We also want to compare to the theoretical prediction:
fp0=data(1,2); %This is f'(0), the dimensionless scaled centerline velocity
centerline_velocity_ratio=ucl/Uc/xmed.^.2/fp0 %Normalized to theoretical value
figure(6)
yplot=[-ylim:.01:ylim]';
plot(yadj(ikeep),uadj(ikeep),'o',yplot,[ones(size(yplot)),yplot.^2]*uf1,...
yplot,[ones(size(yplot)),abs(yplot)]*uf2,0,ucl,'*')
grid on
xlabel('y (cm) adjusted to xmed')
ylabel('u (cm/s) adjusted to xmed')
legend('data','parabolic fit','triangle fit','centerline velocity estimate')
title(['Plot for estimating centerline velocity at x = ',num2str(xmed)])

% Now we estimate the half-width of the profile.  What we shall do is try
% to determine the value of y where the profile is half the centerline
% velocity.  To get a better value, we shall take data points around this
% value, fit a line to the data, and then do some interpolation on the
% line.  Note that this requires that we get the centerline velocity right
% - so that has to be good in order for this procedure to work!

ikeep=find(abs(xc-xmed)<xd&uadj>0.35*ucl&uadj<0.65*ucl);
amat=[abs(yadj(ikeep)),ones(size(yadj(ikeep)))];
hwfit=amat\uadj(ikeep);
halfwidth=(0.5*ucl-hwfit(2))/hwfit(1)
% We also want to compare to the theoretical prediction, which uses the
% results of heated wire.
itheory=find(data(:,2)>0.475*fp0&data(:,2)<0.525*fp0);
amat=[data(itheory,1),ones(size(data(itheory,1)))];
xth=amat\data(itheory,2);
halfwidththeory=(0.5*fp0-xth(2))/xth(1);

halfwidth_ratio=halfwidth/xmed^.4/Ly/halfwidththeory

figure(7)
yplot=[0 halfwidth*2]';
plot(abs(yadj(ikeep)),uadj(ikeep),'o',yplot,[yplot,ones(size(yplot))]*hwfit,'g')
xlabel('|y| (cm) adjusted to xmed')
ylabel('u (cm/s) adjusted to xmed')
title('Plot of data near half max for determining profile half-width')


%Finally, we add the estimated width and magnitude to the dimensional plot:
figure(2)
hold on
plot([halfwidth,-halfwidth],ucl/2*[1,1],'rx',0,ucl,'g*')
hold off
grid on
text(1,ucl/2,char(['centerline velocity = ',num2str(ucl)],['half width = ',num2str(halfwidth)]))
legend(char(labels2,['x = ',num2str(xmed+xd)],['x = ',num2str(xmed-xd)],'scaled theory','half-width est','centerline est'))



