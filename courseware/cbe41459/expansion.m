function texpout=expansion(varargin)
%This function takes in glycerine concentration (mass fraction) and the operating temperature
%and returns the thermal expansion coefficient of the aqueous solution.  The data used
%comes from the Dow website.  The average deviation is 1% from measured values.  Output is in 
%units of 1/C.
%
%The concentration and temperature may be arrays.  If only one variable
%is used, it is taken to be concentration and the results are computed for
%22 degrees C.


%we have the thermal expansion coefficients in 1/C
texp=[0.000615	0.000615	0.00061
0.00062	0.000615	0.000605
0.000615	0.000615	0.000615
0.00061	0.000615	0.00062
0.00062	0.000615	0.00061
0.00058	0.00057	0.000565
0.00054	0.000545	0.00055
0.000485	0.000495	0.00051
0.00043	0.000435	0.000445
0.00037	0.000385	0.0004
0.0003	0.000315	0.000325
0.00023	0.000255	0.00028
0.00018	0.000205	0.00023];

%the concentrations (mass fraction) are
conc=[100
97.5
95
90
80
70
60
50
40
30
20
10
0]/100;

%the given temperatures (C)
temp=[17.5	20	22.5];

%the absolute temperatures
tabs=temp+273.15;

%we do a cubic fit of the thermal expansion coefficient on the concentration
cmodel=[ones(size(conc)),conc,conc.^2,conc.^3];
xcall=cmodel\(texp);

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

%now we obtain the output
cmout=[ones(size(c)),c,c.^2,c.^3];
tmout=[ones(size(t)),(t+273.15)];
texpout=cmout*(tmout*xtall)';

%Notes: The Dow website data is located at the url:
%
%http://www.dow.com/glycerine/resources/table14.htm
%
%This data may be graphically compared to the results
%of the fitting algorithm by copying the data, concentration,
%and temperature from the top of the function, pasting
%it into the command line of Matlab, and then executing
%the following command:
%
%figure;plot(conc,texp,'*',conc,expansion(conc,temp));grid on;zoom on;
%
%D.T.Leighton, University of Notre Dame
%8-26-2010
