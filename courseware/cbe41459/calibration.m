echo on
%this program will calibrate the image.  It takes inputs
%from the mouse and calculates the pixel to centimeter ratio
%which will then be used in pivplot.m
%
%IMPORTANT: do not include spaces in your image's file name!
%
echo off

%first load the image for calibration
calname=input('Enter calibration file name including extension: ','s');
figure(1)
imshow(calname)
title('Calibration Image')

%calculate the pixel distance in the vertical direction
xdist=input('Enter the vertical distance to be calibrated in cm: ');
disp('Click on two points that are separated by the calibration length.');
[y1,x1]=ginput(2); %waits for two points on the image to be clicked
xcalib=abs(x1(1)-x1(2))/xdist  %pixels/cm

%calculate the pixel distance in the horizontal direction
ydist=input('Enter the horizontal distance to be calibrated in cm: ');
disp('Click on two points that are separated by the calibration length.');
[y2,x2]=ginput(2);
ycalib=abs(y2(1)-y2(2))/ydist  %pixels/cm

%calculate the reference height of the center of the image:
xref=input('Enter the vertical position of a reference point in cm: ');
disp('Click on the vertical reference point.');
[y2,x2]=ginput(1);

%We need to check whether the image is widescreen or standard.
[xn yn zn]=size(double(imread(calname)));

xcenter=xref-(xn/2-x2)/xcalib; 

echo on
%Calibration Complete!
%The vertical position of the image center is:

xcenter
%You need to subtract off the vertical position of the wire from this value,
%and then enter it as the vertical distance of the image center from the wire
%in pivplot.m.  You will also need to edit pivplot.m to use the correct x and
%y pixel calibrations!  You are now set up to run pivanal and/or pivplot.
%
echo off


