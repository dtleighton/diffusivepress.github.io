%This script calculates and plots up Uc and Ly at 1cm above the wire for a
%range of concentrations.  We base it on the volume fraction of
%glycerin, then calculate the true mass fraction of glycerin in the
%solution, and then figure out how these variables change.  We use a base
%power of 8 Watts, a typical value, and the default temperature of 22 deg C.

%OK, we start with volume fraction:

vf=[0:.05:.67]'; %modify this to get your composition!

bf=0.997; %a guess as to the true wt% of the glycerin (USP grade)

conc=vf*bf*density(bf)./(vf*density(bf)+(1-vf)*density(0));

q=3.44; %watts input - (change to the value for your experiment!)
leng=11.4; %length of wire in cm - (change to the value for your experiment!)
%OK, we need the kinematic viscosity and thermal diffusivity:

nu=viscosity(conc)./density(conc);

alpha=conductivity(conc)./density(conc)./htcapacity(conc);

%And we get Uc:

uc=(((q/leng)*expansion(conc)*980.0./density(conc)./htcapacity(conc)).^2.0./nu).^0.2;

%and Ly:

ly=(nu.^3.0.*density(conc).*htcapacity(conc)/(q/leng)./expansion(conc)/980).^.2;

tc = uc.*nu./ly.^2.0./expansion(conc)/980;

figure(1)
plot(vf,uc,vf,ly,vf,tc)
xlabel('volume fraction glycerin')
ylabel('Uc (cm/s), Ly (cm) & Tc (C)')
legend('Uc','Ly','Tc')
title(['Velocity magnitude, width, and temperature at ',num2str(q),'W power'])
grid on
zoom on

