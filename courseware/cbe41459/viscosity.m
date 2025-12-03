function viscout=viscosity(varargin)
%This function takes in the concentration (mass fraction) of glycerine and
%the operating temperature and returns the viscosity of the aqueous solution.
%It uses the data from the Dow website, and captures this data to about 1.2%
%average deviation.  The output is in poise.
%
%The concentration and temperature may be arrays.  If only one variable
%is used, it is taken to be concentration and the results are computed for
%22 degrees C.

%We have the Dow data:
temp=[0 10 20   30      40];
conc=[0
10
20
30
40
50
60
65
67
70
75
80
85
90
95
100]/100;
visc=[1.792     1.308   1.005   0.8007  0.656
2.44    1.74    1.31    1.03    0.826
3.44    2.41    1.76    1.35    1.07
5.14    3.49    2.5     1.87    1.46
8.25    5.37    3.72    2.72    2.07
14.6    9.01    6       4.21    3.1
29.9    17.4    10.8    7.19    5.08
45.7    25.3    15.2    9.85    6.8
55.5    29.9    17.7    11.3    7.73
76      38.8    22.5    14.1    9.4
132     65.2    35.5    21.2    13.6
255     116     60.1    33.9    20.8
540     223     109     58      33.5
1310    498     219     109     60
3690    1270    523     237     121
12070   3900    1410    612     284]/100;
tabs=temp+273.16;

%We do a fourth order fit to the log of the viscosity on concentration:
cmodel=[ones(size(conc)),conc,conc.^2,conc.^3,conc.^4];
xcall=cmodel\log(visc);

%Now we do a quadratic fit vs. 1/temperature (abs)

tmodel=[ones(size(tabs')),1.0./tabs',1.0./tabs'.^2];
xtall=tmodel\xcall';

%We identify the target concentrations and temperatures:
c=varargin{1};
if nargin<2;t=22;else;t=varargin{2};end

%We ensure that c and t are column vectors:
[n m]=size(c);
if m>n;c=c';end
[n m]=size(t);
if m>n;t=t';end

%Now we are ready to get the output:
cmout=[ones(size(c)),c,c.^2,c.^3,c.^4];
tmout=[ones(size(t)),1.0./(t+273.16),1.0./(t+273.16).^2];
viscout=exp(cmout*(tmout*xtall)');

%Notes: The Dow website data is located at the url:
%
%http://www.dow.com/glycerine/resources/table18.htm
%
%This data may be graphically compared to the results
%of the fitting algorithm by copying the data, concentration,
%and temperature from the top of the function, pasting
%it into the command line of Matlab, and then executing
%the following command:
%
%figure;semilogy(conc,visc,'*',conc,viscosity(conc,temp));grid on;zoom on;
%
%D.T.Leighton, University of Notre Dame
%8-26-2010
