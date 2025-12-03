function dtout=temprise(varargin)
%This function takes the data on temperature rise of mixing from the dow
%website graph (digitized), fits it to a cubic, and then outputs the temperature
%rise for selected mass fractions of glycerin.  The input may be an array, and 
%the output will be a column vector.  Dow data is for pure glycerin dilution
%at 0 degrees C.  Units are in degrees C.

data =[0         0
    0.1053    1.5000
    0.2035    2.5775
    0.3053    3.5282
    0.4070    4.2887
    0.5018    4.8169
    0.5649    5.1127
    0.6105    4.9648
    0.7123    4.4789
    0.8105    3.7183
    0.9018    2.5352
    1.0000         0];

conc=data(:,1);
deltat=data(:,2);
a=[conc.*(1-conc),conc.^2.0.*(1-conc)];
x=a\deltat;

c=varargin{1};
[n m]=size(c);if n<m; c=c';end; %convert to column vector if necessary
dtout=[c.*(1-c),c.^2.0.*(1-c)]*x;

%Notes: The Dow website data is located at the url:
%
%http://www.dow.com/glycerine/resources/temprise.htm
%
%The graph was digitized manually in matlab to generate the data table.
%
%This data may be graphically compared to the results
%of the fitting algorithm by copying the data, concentration,
%and temperature from the top of the function, pasting
%it into the command line of Matlab, and then executing
%the following command:
%
%figure;plot(conc,deltat,'*',conc,temprise(conc,temp));grid on;zoom on;
%
%D.T.Leighton, University of Notre Dame
%8-26-2010
