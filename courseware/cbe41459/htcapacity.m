function hout=htcapacity(varargin)
%This function takes in the glycerine concentration (mass fraction) and
%temperature and returns the heat capacity of the solution.  The values
%at 15 and 32 degrees C from Perry's (3rd edition) are used.  The output
%is a column vector (or array) in units of J/g/degK.
%
%The concentration and temperature may be arrays.  If only one variable
%is used, it is taken to be concentration and the results are computed for
%22 degrees C.

%heat capacity of glycerine and water (J/g/K)
hcap=[0.99976 0.99797
0.961 0.960
0.929 0.924
0.851 0.841
0.765 0.758
0.67 0.672
0.555 0.576]*4.184;

conc=[0 .1 .2 .4 .6 .8 1]';

%temperature in C
temp=[15 32];
tabs=temp+273.16;

%we do a quadratic fit of the heat capacity on the concentration
cmodel=[ones(size(conc)),conc,conc.^2];
xcall=cmodel\hcap;

%Now we do a linear fit vs. temperature (abs)
tmodel=[ones(size(tabs')),tabs'];
xtall=tmodel\xcall';

%We identify the target concentrations and temperatures:
c=varargin{1};
if nargin<2;t=22;else;t=varargin{2};end

%We ensure that c and t are column vectors:
[n m]=size(c);
if m>n;c=c';end
[n m]=size(t);
if m>n;t=t';end

%We have the output values:
cmout=[ones(size(c)),c,c.^2];
tmout=[ones(size(t)),(t+273.15)];
hout=cmout*(tmout*xtall)';

%Notes: The data is taken from p. 235 of Perry's Handbook, 3rd ed.
%
%The value for pure glycerin is a bit below the range listed on the
%Dow website.
%
%This data may be graphically compared to the results
%of the fitting algorithm by copying the data, concentration,
%and temperature from the top of the function, pasting
%it into the command line of Matlab, and then executing
%the following command:
%
%figure;plot(conc,hcap,'*',conc,htcapacity(conc,temp));grid on;zoom on;
%
%D.T.Leighton, University of Notre Dame
%8-26-2010
