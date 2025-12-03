echo on
%This program takes in the particle identifications and linked lists made by
%the program pivanal.m and generates a quiver plot of velocities.  You need
%to enter the basename of the picture images, as the required data is produced
%by the program pivanal.m, and stored under the name 'basenamelinks.mat'.  The
%velocity data generated here is stored under the name 'basenamevel.mat'.

%You will need to enter in the vertical location of the center of the video images
%(in cm, relative to the wire), as well as the pixel calibration factors 
%in both x (up) and y (across) directions.

%First we load in the data:
echo off

basename=input('Enter the base filename: ','s');

filename=[basename,'links.mat'];
load(filename)

%We now have the three data variables nlinks xlinks and ylinks.  First we
%adjust the x and y center and calibration.  Remember that x should be "up",
%and that the data is in row-column format.

%We have the center of the image and calibration factors (these will have to
%be measured for every camera setting!):

xcalib=257.3; %pixels/cm
ycalib=256.8; %pixels/cm

xcenter=input('Enter the vertical center of the image relative to the wire: ');

%We need to know the size of the pictures used in obtaining the data (whether
%widescreen or standard).  Since the images go pretty much all the way across
%in the y direction, we can test the maximum pixel value for y.

ytest=max(max(ylinks));

nx=480;
ny=640; %default standard values.

if ytest>641
 nx=540;
 ny=960;
end

if ytest>961
 nx=720;
 ny=1280;
end

if ytest>1281
 nx=1080;
 ny=1920;
end

xlinks=((nx-10)/2-xlinks)/xcalib+xcenter;
ylinks=(ylinks-(ny-10)/2)/ycalib;

u=zeros(size(nlinks));
v=zeros(size(nlinks));
xc=zeros(size(nlinks));
yc=zeros(size(nlinks));

for i=1:length(nlinks)
 x=xlinks(i,1:nlinks(i))';
 y=ylinks(i,1:nlinks(i))';
 t=[0:nlinks(i)-1]'/60; %The time base in seconds
 a=[ones(nlinks(i),1),t];
 fitu=a\x;
 fitv=a\y;
 u(i)=fitu(2);
 v(i)=fitv(2);
 xc(i)=mean(x);
 yc(i)=mean(y);
end

% We need to determine the center of the velocity profile in the horizontal
% direction.  This is fairly tricky to actually measure (particularly if
% the profile is unstable!), so instead we calculate it from the velocity
% profile itself.  We divide the data up into strips at different vertical
% positions and then determine the maximum velocity in each strip.  We then
% plot up the lateral location of these max velocity points vs. height and
% fit a line to the data.  This will allow us to automatically correct for
% lateral drift of the velocity profile, but only works if there are a
% sufficient number of data points that we get a good read on the maximum
% in each strip.  Note that you may need to play with the width and the
% number of strips to get better measurement of the center of the profile!

width = 1.5; %The width of the measurements above x=0
nstrips = 6; %The number of strips we look at
loweredge = [0:nstrips-1]/nstrips*width+median(xc)-width/2;
upperedge = loweredge+width/nstrips;
xmax = zeros(nstrips,1);
ymax = zeros(nstrips,1);

for k=1:nstrips
    i = find(xc>loweredge(k)&xc<upperedge(k));
    [junk,ki] = max(u(i));
    xt = xc(i);
    yt = yc(i);
    ymax(k) = yt(ki);
    xmax(k) = xt(ki);
end

%We fit a line to the data:
a = [ones(size(xmax)),xmax];
xfit = a\ymax % The second coefficient is the slope of the drift

% We also need to determine the error in the calculated center of the image
% to make sure that we are getting it right!
k=(a'*a)^-1*a';
sigy = norm(a*xfit-ymax)/(nstrips-2)^.5;
varxfit=k*k'*sigy^2;
sigxfit = diag(varxfit).^.5 % The standard deviation of the fitting coef.

y0 = [1,median(xc)]*xfit %The location of the center at the median
sigy0 = ([1,median(xc)]*varxfit*[1,median(xc)]')^.5 %The standard deviation

figure(5)
plot(xmax,ymax,'o',xmax,a*xfit)
title('Determination of Profile Center')
xlabel('vertical position')
ylabel('horizontal location of maximum')
legend('data',['xfit = ', num2str(xfit')])

% If the uncertainty in y0 is too large, we determine the center by hand.
% This will not allow for a slope to the profile (it is assumed to be
% zero), but it will allow us to get the center "eyeballed" correctly.
if sigy0>0.05 %We choose a value of 0.5mm for the maximum uncertainty.
 figure(5)
 plot(yc,u,'*')
 title('The profile center could not be determined automatically.  Click on the center in the graph.')
 [y0,junk]=ginput(1);
 xfit(1) = y0;
 xfit(2) = 0;
end

%And we subtract this shift off from the original data:
yc = yc - [ones(size(xc)),xc]*xfit;
[nn nm] = size(ylinks);
ylinks=ylinks - ([ones(size(xc)),xc]*xfit)*ones(1,nm);

%There are errors associated with PIV, particularly if there is a spurious image
%in the camera which simply doesn't move.  This yields a zero velocity, which can
%lead to errors.  We can reduce this problem by only keeping data points for which
%the velocity is greater than some threshold fraction of the average velocity:
speed=(u.^2+v.^2).^.5;
minspeed=.1*mean(abs(speed));
i=find(speed>minspeed);
u=u(i);
v=v(i);
yc=yc(i);
xc=xc(i);
xlinks=xlinks(i,:);
ylinks=ylinks(i,:);

%Now we can plot all these calculated velocities up.  We reverse the shift
%so that the quiver plot is in the original image reference frame.
figure(4)
plot(ylinks+([ones(size(xc)),xc]*xfit)*ones(1,nm),xlinks,'.')
title(['Quiver plot of measured velocities for ',basename])
xlabel('y (cm)')
ylabel('x (cm)')
hold on
quiver(yc+[ones(size(xc)),xc]*xfit,xc,v,u,'o')
hold off
zoom on

%We can also plot up a velocity profile from this data.  Let's just use data
%close to the median of the vertical location:
medx = median(xc);
i=find(xc>(medx-1)&xc<(medx+1)); %We look +- 1.0 cm around the median
figure(6)
plot(yc(i),u(i),'*')
xlabel('y (cm)')
ylabel('u (cm/s)')
title(['Velocity Profile for ',basename,', x = ',num2str(medx),'+-1.0 cm, xfit = ', num2str(xfit')])
grid

%The y-velocity is also of interest.
[ysort,k]=sort(yc(i));
vsort=v(i(k));
figure(7)
plot(ysort,vsort,'*')
title(['Horizontal Velocity Profile for ',basename,', x = ',num2str(medx),'+-1.0 cm, xfit = ', num2str(xfit')])
xlabel('y (cm)')
ylabel('v (cm/s)')
grid

%Finally, we save the velocity data and the centroids of the links as a data file
%for later analysis:
outname=[basename,'vel.mat'];
save(outname,'xc','yc','u','v','xfit')




