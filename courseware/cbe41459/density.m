function denout=density(varargin)
%This function takes in glycerine concentration (mass fraction) and the operating temperature
%and returns the density of the aqueous solution.  The data used comes from the Dow website.  
%The maximum deviation between the measured and fitted values is 0.1%.
%
%The concentration and temperature may be arrays.  If only one variable
%is used, it is taken to be concentration and the results are computed for
%22 degrees C.

%we have the density in g/cm^3 from the DOW website.
den=[1.26415	1.26108	1.25802	1.25495
1.2381	1.2351	1.232	1.2289
1.2116	1.2085	1.20545	1.2024
1.18415	1.18125	1.1784	1.17565
1.1565	1.1538	1.15105	1.1483
1.1287	1.1263	1.12375	1.1211
1.10145	1.0993	1.0971	1.09475
1.07455	1.0727	1.0707	1.06855
1.0484	1.0469	1.04525	1.0435
1.02325	1.0221	1.0207	1.01905
0.99913	0.99823	0.99708	0.99568];

%the concentrations in mass fractions are
conc=[100
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

%the range of temperatures (C) is as follows:
temp=[15	20	25	30];

%the temperature in Kelvin
tabs=temp+273.15;

%we do a quadratic fit of the reciprocal of the density on the concentration
cmodel=[ones(size(conc)),conc,conc.^2];
xcall=cmodel\(1.0./den);

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
cmout=[ones(size(c)),c,c.^2];
tmout=[ones(size(t)),(t+273.15)];
denout=1.0./(cmout*(tmout*xtall)');

%Notes: The Dow website data is located at the url:
%
%http://www.dow.com/webapps/lit/litorder.asp?filepath=glycerine/pdfs/noreg/115-00656.pdf&pdf=true
%
%This data may be graphically compared to the results
%of the fitting algorithm by copying the data, concentration,
%and temperature from the top of the function, pasting
%it into the command line of Matlab, and then executing
%the following command:
%
%figure;plot(conc,den,'*',conc,density(conc,temp));grid on;zoom on;
%
%D.T.Leighton, University of Notre Dame
%8-26-2010
