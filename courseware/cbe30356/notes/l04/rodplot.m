%% Temperature Profile of a Cooled Rod
% We have the dimensionless temperature profile, which is a function of
% lambda.  We are interested in the profile for values of lambda ranging
% from 1 to 3.

z = [0:.01:1];

la = [1,1.35,2.2,3]';

T = cosh(la*z) - ((sinh(la)./cosh(la))*ones(size(z))).*sinh(la*z);

figure(1)
plot(z,T)
grid on
xlabel('z')
ylabel('T^*')
title('Temperature Profile for Different Values of Lambda')
legend(str2mat(num2str(la)))