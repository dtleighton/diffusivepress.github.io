function thermout=conductivity(varargin)
%This function takes in glycerin concentration (mass fraction) and the operating temperature
%and returns the thermal conductivity of the aqueous solution.  The data used
%comes from the Dow website.  The average deviation of the fitted conductivities from
%the measured values is less than 0.4%.
%
%Conductivities are output as an array in units of J/(cm*s*C).
%
%The concentration and temperature may be arrays.  If only one variable
%is used, it is taken to be concentration and the results are computed for
%22 degrees C.

%We have the Dow thermal conductivity data in cal/(cm*s*C)
therm=[0.00138	0.00141	0.00145	0.00149
0.00133	0.00137	0.0014	0.00144
0.0013	0.00133	0.00137	0.0014
0.00125	0.00128	0.00131	0.00134
0.00121	0.00124	0.00127	0.00129
0.00117	0.00119	0.00122	0.00125
0.00112	0.00115	0.00117	0.0012
0.00109	0.00111	0.00114	0.00116
0.00105	0.00107	0.00108	0.0011
0.00102	0.00103	0.00105	0.00106
0.00097	0.00099	0.001	0.00101
0.00094	0.00095	0.00096	0.00098
0.0009	0.00091	0.00091	0.00092
0.00086	0.00087	0.00088	0.00089
0.00084	0.00084	0.00085	0.00085
0.0008	0.00081	0.00081	0.00081
0.00077	0.00078	0.00078	0.00078
0.00074	0.00074	0.00074	0.00074
0.00072	0.00072	0.00072	0.00072
0.0007	0.0007	0.0007	0.0007
0.00068	0.00068	0.00068	0.00068];

%we convert to Joules/(cm*s*C):
therm=therm*4.184;

%the temperatures in degrees celsius
temp=[10	20	30	40];

%the measured glycerine concentrations
conc=[0
5
10
15
20
25
30
35
40
45
50
55
60
65
70
75
80
85
90
95
100]/100;

%absolute temperature
tabs=temp+273.15;

%we do a quadratic fit of the thermal conductivity on the concentration
cmodel=[ones(size(conc)),conc,conc.^2,conc.^3];
xcall=cmodel\therm;

%Now we do a linear fit vs. 1/temperature (abs)
tmodel=[ones(size(tabs')),1.0./tabs'];
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
tmout=[ones(size(t)),1.0./(t+273.15)];
thermout=cmout*(tmout*xtall)';

%Notes: The Dow website data is located at the url:
%
%http://www.dow.com/glycerine/resources/table13.htm
%
%This data may be graphically compared to the results
%of the fitting algorithm by copying the data, concentration,
%and temperature from the top of the function, pasting
%it into the command line of Matlab, and then executing
%the following command:
%
%figure;plot(conc,therm,'*',conc,conductivity(conc,temp));grid on;zoom on;
%
%D.T.Leighton, University of Notre Dame
%8-26-2010

