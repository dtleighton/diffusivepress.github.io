%% Lecture 1 Problem of the Day Solution
% In this script we calculate the thermal conductivity of N2 and compare it
% to literature values.  We use the Eucken formula for low density gases,
% which relates the viscosity (calculated from the Chapman Enskog model) to
% the thermal conductivity via the heat capacity.

%% Heat Capacity Measurements
% The thermal conductivity requires knowledge of a number of molecular
% parameters and also the heat capacity as a function of temperature.  This
% last is needed because a diatomic (or polyatomic) gas has additional
% thermal modes which enable them to carry energy other than just pure
% kinetic energy, and these become active at different temperatures.  Thus,
% the heat capacity goes up with temperature.  From the data, it is pretty
% constant up to around 400°K where vibrational modes start to become
% active.  Let's look at this:

% From the engineering toolbox cut and paste:
table1=[175	1.039
200	1.039
225	1.039
250	1.039
275	1.039
300	1.040
325	1.040
350	1.041
375	1.042
400	1.044
450	1.049
500	1.056
550	1.065
600	1.075
650	1.086
700	1.098
750	1.110
800	1.122
850	1.134
900	1.146
950	1.157
1000	1.167
1050	1.177
1100	1.187
1150	1.196
1200	1.204
1250	1.212
1300	1.219
1350	1.226
1400	1.232
1500	1.244
1600	1.254
1700	1.263
1800	1.271
1900	1.278
2000	1.284
2100	1.290
2200	1.295
2300	1.300
2400	1.304
2500	1.307
2600	1.311
2700	1.314
2800	1.317
2900	1.320
3000	1.323
3500	1.333
4000	1.342
4500	1.349
5000	1.355
5500	1.362
6000	1.369];
t = table1(:,1);
cp = table1(:,2);

figure(1)
semilogx(t,cp,'o-')
grid on
xlabel('T, °K')
ylabel('specific heat capacity, kJ/(kg °K)')
title('Specific Heat Capacity at 1 bar for N2')

%% Calculation of Viscosity
% The Chapman-Enskog model for the viscosity needs the molecular weight,
% Lennard-Jones parameters and collision integrals.  We can get these from
% the tables E-1 and E-2 from BSL:

sigma = 3.667e-8; %units of cm
eok = 99.8; %units of °K

% and the formula:
omega = @(t) 1.16145./(t/eok).^0.14874 + 0.52487./exp(0.77320*t/eok)+2.16178./exp(2.43787*t/eok)

% The other parameters:
m = 28.013/6.02214e23 %g/molecule
kb = 1.380649e-16 %units in cgs

% and the Chapman-Enskog equation:
visc = @(t) 5/16*(pi*m*kb*t).^.5./(pi*sigma^2*omega(t))

% Now we compare to the data from the engineering toolbox to make sure
% we've done it right:

table2=[80	-193.2	-315.7	0.1	1	5.623	0.00562	0.1174	3.778	0.01360	1.284	13.82
100	-173.2	-279.7	0.1	1	6.958	0.00696	0.1453	4.676	0.01683	2.024	21.79
120	-153.2	-243.7	0.1	1	8.244	0.00824	0.1722	5.540	0.01994	2.903	31.25
140	-133.2	-207.7	0.1	1	9.480	0.00948	0.1980	6.370	0.02293	3.911	42.10
160	-113.2	-171.7	0.1	1	10.67	0.01067	0.2228	7.170	0.02581	5.043	54.28
180	-93.2	-135.7	0.1	1	11.81	0.01181	0.2467	7.936	0.02857	6.289	67.69
200	-73.2	-99.7	0.1	1	12.91	0.01291	0.2696	8.675	0.03123	7.648	82.32
220	-53.2	-63.7	0.1	1	13.97	0.01397	0.2918	9.387	0.03379	9.107	98.03
240	-33.2	-27.7	0.1	1	15.00	0.01500	0.3133	10.08	0.03629	10.68	114.9
260	-13.2	8.3	0.1	1	15.99	0.01599	0.3340	10.74	0.03868	12.33	132.7
280	6.9	44.3	0.1	1	16.96	0.01696	0.3542	11.40	0.04103	14.09	151.6
300	26.9	80.3	0.1	1	17.89	0.01789	0.3736	12.02	0.04328	15.93	171.5
320	46.9	116.3	0.1	1	18.80	0.01880	0.3926	12.63	0.04548	17.85	192.2
340	66.9	152.3	0.1	1	19.68	0.01968	0.4110	13.22	0.04761	19.86	213.8
360	86.9	188.3	0.1	1	20.55	0.02055	0.4292	13.81	0.04971	21.96	236.4
400	126.9	260.3	0.1	1	22.21	0.02221	0.4639	14.92	0.05373	26.37	283.9
500	226.9	440.3	0.1	1	26.06	0.02606	0.5443	17.51	0.06304	38.69	416.4
600	326.9	620.3	0.1	1	29.58	0.02958	0.6178	19.88	0.07156	52.70	567.2
700	426.9	800.3	0.1	1	32.83	0.03283	0.6857	22.06	0.07942	68.24	734.5
800	526.9	980.3	0.1	1	35.89	0.03589	0.7496	24.12	0.08682	85.25	917.6
900	626.9	1160.3	0.1	1	38.78	0.03878	0.8099	26.06	0.09381	103.63	1115.5
1000	726.9	1340.3	0.1	1	41.54	0.04154	0.8676	27.91	0.1005	123.34	1327.6];

% We extract the relevant columns:
tvisc = table2(:,1);
viscdata = table2(:,7)/100; %we convert to poise (cgs)

figure(2)
semilogx(tvisc,viscdata,'o-',tvisc,visc(tvisc),'--')
axis([80 1000 0 5e-4])
grid on
xlabel('T, °K')
ylabel('viscosity (poise)')
legend('engineering toolbox data','Chapman-Enskog model')
title('Comparison of the viscosity of N2 with the Chapman-Enskog model')

%% Calculation of the Thermal Conductivity
% Now we put these two together to get the thermal conductivity.  We need
% the additional parameter R/M from the molecular weight and the ideal gas
% law constant.  We also have to be careful of units: from the table Cp is
% in KJ/kg°K so we would need to multiply by 1000 to get SI, and the
% viscosity is in cgs (poise) so we divide by 10 to get it into SI (Pa*s).
% We can only evaluate the Eucken Formula where we have heat capacity data.

rom = 8.3145/28.013; % R/M in the units of kJ/(kg °K)

condeucken = visc(t).*(cp+5/4*rom)/10*1000; % units of W/(m°K)

% and we have the data from the engineering toolbox again:

table3=[78	-195	7.253	0.00624	78	-320	0.00416	0.00619
98	-175	9.216	0.00792	89	-300	0.00479	0.00713
123	-150	11.58	0.00995	116	-250	0.00633	0.00942
148	-125	13.84	0.01190	144	-200	0.00780	0.01161
173	-100	16.02	0.01378	172	-150	0.00920	0.01370
198	-75	18.12	0.01558	200	-100	0.01055	0.01570
223	-50	20.15	0.01732	228	-50	0.01184	0.01763
231	-42	20.78	0.01787	231	-44	0.01200	0.01785
248	-25	22.10	0.01901	244	-20	0.01260	0.01875
263	-10	23.25	0.01999	255	0	0.01309	0.01948
268	-5	23.62	0.02031	266	20	0.01358	0.02021
278	5	24.37	0.02095	278	40	0.01406	0.02092
283	10	24.74	0.02127	283	50	0.01429	0.02127
293	20	25.47	0.02190	294	70	0.01476	0.02197
298	25	25.83	0.02221	300	80	0.01500	0.02232
303	30	26.19	0.02252	305	90	0.01523	0.02266
323	50	27.62	0.02375	311	100	0.01546	0.02300
348	75	29.35	0.02524	339	150	0.01658	0.02468
373	100	31.04	0.02669	366	200	0.01768	0.02631
398	125	32.69	0.02811	394	250	0.01874	0.02789
423	150	34.30	0.02949	422	300	0.01978	0.02943
448	175	35.87	0.03085	478	400	0.02178	0.03241
473	200	37.42	0.03217	533	500	0.02369	0.03526
573	300	43.32	0.03725	644	700	0.02733	0.04066
623	350	46.13	0.03966	700	800	0.02906	0.04325
673	400	48.86	0.04202	755	900	0.03075	0.04576
773	500	54.14	0.04655	811	1000	0.03240	0.04821];
tcond = table3(:,1);
conddata = table3(:,3)/1000; % converting units to SI

% and we plot it up:
figure(3)
semilogx(tcond,conddata,'o-',t,condeucken,'--')
grid on
xlabel('T, °K')
ylabel('Thermal Conductivity (W/m°K)')
title('Thermal Conductivity of N2 at 1 bar')
legend('data','Eucken Formula')
axis([80 1000 0 0.06])

%% Conclusion
% Examination of the graphs shows that the Chapman-Enskog model matches the
% Engineering Toolbox data pretty well for the viscosity over a broad range
% of temperatures, and the Eucken formula gets the thermal conductivity
% within 5% for temperatures below 300°K.  The predictions deviate a bit at
% higher temperature, but are still within 10% at 770°K, the highest
% temperature for which data is given.  Looking at the graph, we see that k
% is a weak function of temperature, and examining the other graphs and
% data from the engineering toolbox we see that it is pretty much
% independent of pressure until you get up to 100 bar or more, also as
% expected.  The low density gas formula breaks down as you approach the
% critical point, which for N2 is 34 bar and 126°K.