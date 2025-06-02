%% Problem of the Day 06: Heating of a Slab from Both Sides
% We plot up the dimensionless profile of the temperature distribution for
% a slab with heating at both sides.  We take t to be a scalar and x to be
% a row vector.  That way the product of a column vector n and row vector x
% produces a matrix which can be summed down the columns (the default for
% matlab) to sum the series without using a loop.  Because the lead
% eigenvalue is pi, the decaying solution vanishes as exp(-pi^2*t), so it
% is all over by a dimensionless time of around 0.25 or so.

Tinf = @(t,x) t + 0.5*x.^2 - 1/6

n = [1:100]'; %We keep 100 eigenvalues (overkill)
Td = @(t,x) sum(2*(-(-1).^n./(n*pi).^2.*exp(-(n*pi).^2*t))*ones(size(x)).*cos(n*pi*x))

x = [-1:.01:1];

figure(1)
tall = [0,.05,.1,.2,.4]';
colors = 'rgbcmyk';

for i = 1:length(tall)
    plot(x,Tinf(tall(i),x)+Td(tall(i),x),colors(i))
    hold on
end
for i = 1:length(tall)
    plot(x,Tinf(tall(i),x),[colors(i),'--'])
    hold on
end

hold off
xlabel('x')
ylabel('T')
legend(num2str(tall))
grid on
title('Temperature Profile of a Slab (dashed lines are Tinf)')
