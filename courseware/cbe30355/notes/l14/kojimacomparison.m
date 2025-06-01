%% Comparison of Kojima Drop Data with Pignatel Particle Cloud Data
% In this script we compare the miscible torus expansion rate measured by
% Kojima et al. with the suspension drop breakup data of Pignatel et al. in
% an effort to determine conditions under which breakup may be more easily
% visualized in a classroom demonstration.  The droplet breakup phenomenon
% is visually appealing and can be used to motivate students into further
% study of phenomena in fluid mechanics.
%
% The fate of a miscible drop settling through a fluid at low Re is
% complex, but appears to be divided into three phases.  Under creeping
% flow conditions a spherical drop is neutrally stable even in the absence
% of surface tension. It will remain spherical, as the Stokes flow stresses
% on all parts of the surface are uniform.  In practice, however, it will
% deform due to either the convection of any disturbances to the back of
% the sphere, eventually punching through and making a torus, or via
% inertial forces leading to the formation of an oblate spheroid which
% similarly progresses to a torus.
%
% Once a torus has formed, the next stage of the evolution is expansion. If
% the torus is symmetric, then it will not expand under creeping flow
% conditions due to Stokes flow reversibility.  Machu, et al., however,
% demonstrated that a torus composed of a dilute suspension of negatively
% buoyant particles in the same (viscous) fluid would expand if it were
% asymmetric in the flow (settling) direction.  Inertial forces were
% demonstrated to lead to expansion by Kojima et al.  In their work they
% used asymptotics to calculate the rate of expansion, and showed that it
% should be proportional to the Reynolds number.  When compared to
% measurements, however, they found that a torus expands much more slowly
% than predicted. The discrepancy was attributed to a transient surface
% tension between the two miscible liquids.
%
% In the final stages of expansion the torus becomes unstable and falls
% apart into two (or more) drops which then may form tori that become
% unstable in turn.  This is the principal mechanism for drop breakup in
% the absence of surface tension or shear at low Re.
%
% In recent years a number of investigators have looked at drops made up of
% dilute suspensions of particles.  For small particles the drop velocity
% is much greater than the sedimentation velocity of the particles, and
% thus the cloud or blob behaves as if it were a miscible drop.  Such drops
% have been studied both computationally and experimentally.  The extensive
% experimental measurements of Pignatel et al., for dilute drop suspensions
% are particularly useful for comparison to the experimental measurements
% of Kojima et al for miscible drops.

%% Kojima Drop Measurements
% In these experiments drops of Karo syrup and water were deposited into a
% solution of Karo syrup and water at a different concentration.  Thus,
% both density difference and viscosity ratio were varied, as well as drop
% size.  While drops were observed to form a torus, expand, and break up,
% quantitative measurements were only presented for the expansion phase.
% The characteristic time for expansion would be the ratio of the torus
% major radius to the expansion rate db/dt.  Because the paper presents
% both the major radius b as well as the torus circle radius eps*b, we can
% calculate the initial drop volume and Stokes sedimentation time scale
% based on an undeformed drop.  We can also compute the undeformed drop
% Reynolds number Rec based on the Stokes sedimentation velocity.  The
% Reynolds numbers (calculated in this manner) varied from about .2 to 1,
% while the dimensionless expansion time varied from about 100 to 200. This
% non-dimensionalization is chosen so that it aligns with that used by
% Pignatel.  The data is presented in the same order as Table 2 of that
% paper.

%The fluid viscosity
muf=[.51,.51,.61,.61,.61,.61,.77,.77]';

%The viscosity ratio
lambda = [11.6,7.7,7.1,7.1,4.0,4.0,3.9,2.8]';

%The ring radius
b = [.27,.26,.31,.28,.31,.19,.27,.25]';

%The torus aspect ratio
ep = [.37,.42,.41,.39,.44,.49,.35,.40]';

%We calculate the volume of the torus to get a measure of the original drop
%volume.  This assumes that the ring is symmetric and circular.
v = 2*pi^2*ep.^2.0.*b.^3;

%We have the initial drop radius:
R0 = (v/4/pi*3).^(1/3);

%The measured fluid density
rhof = [1.264,1.264,1.274,1.274,1.274,1.274,1.287,1.287]';

%The measured density difference
drho = [.071,.065,.057,.057,.040,.040,.041,.031]';

g = 980;

%We calculate the Stokes sedimentation velocity:
Us = 2/3*drho*g.*R0.^2./muf./((2+3*lambda)./(1+lambda));

%We calculate the Reynolds number:
Rec = Us.*rhof.*R0./muf

%We have the measured velocity
u = [.93,.90,.95,.75,.82,.42,.36,.32]';

%And the calculated velocities of the tori
ucalc = [0.95 1.0 1.02 0.81 0.88 0.47 0.41 0.34]';

%The measured ring expansion rate
dbdt = [.018,.015,.017,.011,.013,.005,.004,.003]';

%The dimensionless ring expansion time:
te = b./dbdt.*Us./R0

%% Sedimentation Velocity
% It is interesting to compare the measured sedimentation velocities of the
% tori to the Stokes sedimentation velocity of a spherical drop of the same
% volume.  Because of the aspect ratio, it would be expected to be smaller,
% as is observed. The difference between the velocities is due to the
% aspect ratio of the ring yielding greater drag than a drop of the same
% volume but of spherical shape as well as a correction due to drop
% inertia. Kojima calculated the expected velocities using the measured
% aspect ratio of the tori and the first inertial correction.  These are
% plotted as well, and closely match the measured velocities.  The bulk of
% the correction is due to the aspect ratio of the tori as the inertial
% correction is small for these conditions.   The ratio of the measured
% velocities to spherical Stokes drop velocities is 0.59 with a sample
% standard deviation of 0.059.
figure(1)
plot(Us,u,'o',[0 max(Us)],[0 max(Us)],Us,ucalc,'x')
xlabel('Stokes sedimentation velocity (cm/s)')
ylabel('Measured sedimentation velocity (cm/s)')
grid on
legend('data','Us','Calculated velocity')

averatio = mean(u./Us)
stdratio = std(u./Us)
%% Ring Expansion Time
% The dimensionless ring expansion time should be a function of Rec, the
% characteristic Reynolds number for an undeformed drop.  This is plotted
% below.  We fit the data to a power law.  Including all the data, the
% expansion time is a decreasing function of Re, with an exponent of -0.23.
% If we exclude the outlier (this corresponded to the smallest ring radius
% and largest value of epsilon - a small, fat torus), then the slope is
% close to -1/3.  This is intriguing, because it would yield a
% characteristic settling height for ring expansion which is independent of
% drop volume.  The error bars are the 30% uncertainty from the paper by
% Kojima.  Note that there may also be a dependence on the viscosity ratio,
% however because both viscosity ratio and density difference vary
% simultaneously, there is no way to test this from the data.  The
% theoretical calculations of Kojima et al., however, suggests that, much
% like the Stokes sedimentation velocity, the rate of expansion is only a
% weak function of the viscosity ratio.

ak = [ones(size(Rec)),log(Rec)];
x = ak\log(te)

figure(2)
loglog(Rec,te,'o',Rec,exp(ak*x),Rec,113./Rec.^(1/3))
hold on
errorbar(Rec,te,te*0.3,'.')
hold off
axis([.18 1.1 70 300])
xlabel('Re_c')
ylabel('t_e*')
grid on
legend('data','113/Re_c^{0.23}','113/Re_c^{1/3}')

%% Comparison with Suspension Drops
% The data from Pignatel et al for dilute suspensions of particles is given
% below.  This data is for volume fractions of 2% to 10%, and where the
% drop is in the "macro-inertial" regime: the triangles in figure 12a of
% that paper.  This is the appropriate data to use to compare to a miscible
% fluid drop, as the other data in figure 12a would be for the
% "micro-inertial" regime for inertia mediated interactions between
% individual particles.
%
% The time presented by Pignatel is the dimensionless breakup time rather
% than the characteristic ring expansion time, however the two data sets
% are very similar.  Suspension drops were observed to break up when the
% aspect ratio approached 3.  The slope of the suspension drop breakup time
% is slightly greater than that observed for miscible drops, yielding a
% power law fit of 98/Re^0.45, however it again is quite close to the 1/3
% power law yielding heights independent of drop volume.  The 1/3 power law
% fit is added to show the comparison.

greendata=[0.0522  612.5290
    0.0190  551.2761
    0.0584  476.1021
    0.0390  466.3573
    0.0278  413.4571
    0.0558  410.6729
    0.0731  395.3596
    0.0653  387.0070
    0.0522  320.1856
    0.0764  311.8329
    0.0957  306.2645];

leftreddata=[0.1238  306.4935
    0.1187  279.2208
    0.1403  254.5455
    0.1922  246.7532
    0.1346  240.2597
    0.1625  206.4935
    0.1731  206.4935
    0.2226  183.1169
    0.2984  164.9351
    0.2922  155.8442
    0.2134  142.8571
    0.2471  132.4675
    0.1495  120.7792
    0.2922  123.3766
    0.3758  125.9740
    0.4444  109.0909
    0.4262  100.0000
    0.3178  101.2987
    0.2370   85.7143
    0.3178   77.9221
    0.2471   49.3506];

toprightreddata=[0.5540  490.1235
    0.4113  279.0123
    0.5220  260.4938
    0.3724  254.3210
    0.4634  223.4568
    0.4032  214.8148
    0.7031  191.3580
    0.9661  176.5432
    0.7031  166.6667];

bottomrightreddata=[1.4266   14.8855
    1.3251   53.8168
    1.4266   54.9618
    1.3498   70.9924
    1.1017   69.8473
    1.9168   75.5725
    1.4803   82.4427
    1.4266   82.4427
    1.1862   81.2977
    0.8993   79.0076
    0.8828   83.5878
    1.0233   91.6031
    1.2537   90.4580
    0.7074   93.8931
    0.8353   96.1832
    1.0424  111.0687
    0.8667  119.0840
    0.7340  114.5038
    0.6103  117.9389
    0.6450  128.2443
    0.6693  131.6794
    0.7477  133.9695
    0.8200  139.6947
    1.0816  136.2595
    0.9331  145.4198
    0.7340  146.5649
    0.6945  145.4198
    0.5169  147.7099
    0.7074  168.3206
    0.9862  177.4809
    0.8200  113.3588];

alldata=[greendata;leftreddata;toprightreddata;bottomrightreddata];

recp = alldata(:,1);
tb = alldata(:,2);

ap=[ones(size(recp)),log(recp)];
xpignatel=ap\log(tb)

figure(3)
loglog(recp,tb,'ob',recp,98./recp.^(.4536),'b',recp,121./recp.^(1/3),'c')
hold on
loglog(Rec,te,'*r')
errorbar(Rec,te,te*0.3,'.')
hold off
xlabel('Re_c')
ylabel('t_b*, t_e*')
grid on
axis([.01 3 30 1000])
legend('Pignatel data, t_b*','98/Re_c^{.4536}','121/Re^{1/3}','Kojima Data, t_e*')

%% Breakup Height Prediction
% Putting all this together, we can make a prediction of the height
% necessary for drop breakup.  If we use the 121/Re^(1/3) fit, then the
% breakup height is independent of drop volume, and only depends on the
% fluid parameters as the drop radius cancels out.  Because a torus settles
% more slowly than an equivalent drop, the actual height the drop would
% fall to breakup is reduced from this value, with the ratio changing as
% the torus expands.  If we take the ratio of 0.59 from the Kojima
% experiments as characteristic, then the expression for breakup height of
% a suspension drop (viscosity ratio of 1) would be given by:
%
% Hb = 111*(muf.^2./(drho*g.*rhof)).^(1/3)
%
% These values for the Kojima experiments are given below (in cm) and range
% from 16 to 28 cm.  The drop expansion and breakup depicted in the
% pictures of figures 1-4 of Kojima, et al. is the second of these heights
% with a value of 16.4 cm. No heights are reported in the paper, however
% the vessel liquid height is reported to be 82 cm and no photographs were
% taken less than 30 cm from the bottom to avoid wall effects, thus the
% height corresponding to the breakup event depicted in figure 4 had to be
% less than 50 cm.  In addition, for at least some of these conditions (not
% stated) secondary cascade breakup was also observed.  Thus, it is likely
% that the heights calculated below are consistent with the Kojima
% experiments.  Breakup heights are not directly reported by Pignatel
% either, however their vessels were 100 cm in height and multiple breakup
% cascades were also observed under some conditions.

Hb = 111*(muf.^2./(drho*g.*rhof)).^(1/3)

%% Conclusion
% From recent publications, it is apparent that many different processes
% control the ring formation and breakup for falling drops, whether of
% miscible fluids or suspensions.  In particular, suspension drops have
% been shown to break up via purely Stokesian interactions (due to particle
% loss and asymmetry), due to "macro-inertial" effects on the length scale
% of the drop, and due to "micro-inertial" effects on the particle length
% scale.  For fluid drops the breakup is attributed to inertia, however
% there is also the possibility of transient surface tension effects due to
% dissimilar materials (necessary to get a density difference).  The close
% agreement between the expansion time of Kojima and the breakup time of
% Pignatel (where there can be no surface tension effects) make this less
% certain, however.  In addition, while other work has demonstrated the
% existence of such a transient surface tension, it would be expected to
% play more of a role for small drops rather than large ones (e.g., smaller
% Re rather than larger Re).  This is not apparent from the data, where to
% get agreement with ring expansion rates surface tensions (chosen as an
% adjustable parameter) had to be decreased by an order of magnitude for
% smaller drops of the same fluid pairs.
%
% It is clear (and expected) that the breakup time would be a decreasing
% function of Re due to the increasing effect of inertia.  Why it would
% have the observed scaling lying between Re^(-.23) and Re^(-.45) is
% uncertain, however it is usefully approximated by an empirical value of
% Re^(-1/3) which yields a breakup height roughly independent of drop
% volume.  This slightly underpredicts the time observed by Pignatel at low
% Re and overpredicts it at higher Re, but falls within the scatter of the
% data and closely matches that of Kojima.  We shall thus use this
% empiricism in determining optimal fluid/drop combinations.  Note that the
% way in which the drop is introduced also likely affects breakup: the
% significant inertia of a drop falling into a fluid affects the initial
% conditions substantially. The original work of Thomson (1885) found that
% the best rings were produced for drops falling from a height of 1 to 3
% inches.  In the experiments of Kojima drops were released from a height
% of 5 cm to yield an oblate spheroid upon impact.  In the case of
% Pignatel, drops were injected directly into the fluid, likely producing a
% substantially different intial condition.
%
% In order to make a clear demonstration of the phenomenon, it is necessary
% to have a reasonably short breakup height.  This is particularly true if
% it is desired to see a cascade of drop breakup events.  Thus, to make
% things work it is necessary to have a reasonably large density difference
% and low fluid viscosity.  Of the two, the dependence on fluid viscosity
% is somewhat greater, however too low a fluid viscosity (or too high a
% density difference) would lead to velocities and Reynolds numbers well
% beyond the conditions explored by Kojima or Pignatel.  For a
% demonstration in a graduated cylinder, the behavior is further
% complicated by wall reflections.  The diameter of a 500ml cylinder, for
% example, is about 4.5cm inside diameter, thus for a 1 cm diameter drop
% the aspect ratio would be less than 1:5 even before expansion into a
% torus.  This would reduce the sedimentation velocity, but also may limit
% ring expansion and instability.  It also introduces another length scale
% into the problem, and would certainly lead to variations of breakup
% height for different volume drops.
%
% A convenient mixture would be glycerin/water solutions of different
% concentrations.  The glycerin concentration of the fluid would need to be
% 75% by mass for a viscosity of 0.28 poise, and the concentration of the
% drop phase would need to be higher, up to about 95% for a viscosity ratio
% of 10.  Calculated breakup heights based on this combination are depicted
% below, yielding heights ranging from 12 to 24 cm.  A fluid composition of
% 75wt% (71.5% by volume) and a droplet composition of 88wt% (85.6% by
% volume) would be a good combination, yielding a viscosity ratio of 4 and
% a density difference of 0.036 g/cm^3, with a predicted breakup height of
% 14 cm.  This is significantly less than the height of about 35cm
% achievable in a 500ml graduated cylinder. In contrast, if a somewhat
% higher glycerin composition for the fluid is used (e.g., a mass fraction
% of 0.88, yielding a viscosity of 1.1 poise), then the predicted breakup
% height would exceed the height of a 500ml cylinder for all drop glycerin
% concentrations.
%
% For this choice of fluid and drop compositions, a 0.1 ml drop would have
% a diameter of 0.58 cm, a spherical drop Stokes velocity of 2.5 cm/s, and
% a Rec of 3.2.  The latter is somewhat greater than the values explored by
% Kojima, and just a bit larger than suspension drops (in the
% macro-inertial range) examined by Pignatel.  It would be expected to
% break up in about 10 seconds, all reasonable values for a demonstration.
% A 500ml graduated cylinder may be of sufficient height to observe a
% secondary breakup cascade as well.  A smaller drop would put it in the
% same Rec range as those of Kojima, who used drops as small as 0.03 ml,
% however too small a drop becomes both difficult to see in a demonstration
% and to administer in a controlled manner without more extensive equipment
% than dripping from a tube. More viscous base fluids would reduce the
% Reynolds number as well, but at the cost of increasing the height
% required for breakup.  A lower glycerin composition drop would have a
% similar effect due to the smaller density difference.

cf = 0.75;
cd = [cf+.02:.01:0.95]';

hpred = 111*(viscosity(cf)^2./((density(cd)-density(cf))*g*density(cf))).^(1/3);
figure(4)
plot(cd,hpred)
xlabel('drop glycerin mass fraction')
ylabel('predicted breakup height (cm)')
title(['Predicted breakup height for fluid glycerin concentration of ',num2str(cf)])
grid on

%% References
%
% Masami Kojima, E. J. Hinch, and Andreas Acrivos, "The formation and
% expansion of a toroidal drop moving in a viscous fluid", Physics of
% Fluids 27, pp. 19-32 (1984).
%
% Florent Pignatel, Maxime Nicolas, Elisabeth Guazzelli, "A falling cloud of
% particles at a small but finite Reynolds number", J. Fluid Mech. 671, pp.
% 34-51 (2011)
%
% Gunther Machu, Walter Meile, Ludwig C. Nitsche and Uwe Schaflinger,
% "Coalescence, torus formation and breakup of sedimenting drops:
% experiments and computer simulations", J. Fluid Mech. 447, pp. 299-336
% (2001).
%
% Thomson, J. J. & Newall, H. F. "On the formation of vortex rings by drops
% falling into liquids, and some allied phenomena", Proc. R. Soc. Lond. 39,
% pp. 417-435 (1885).
