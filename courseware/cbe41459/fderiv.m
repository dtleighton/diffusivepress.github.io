function fdot=fderiv(eta,f)
%This function provides the derivatives of the vector f.  We pass in
%the prandtl number as a global variable.  The position variable eta
%is not used, but is required for the automatic ode23 solver syntax.
global pr

fdot=zeros(size(f));
fdot(1)=f(2);
fdot(2)=f(3);
fdot(3)=-f(4)+0.2*f(2).^2-0.6*f(1).*f(3); %The x-momentum equation
fdot(4)=-0.6*pr*f(1).*f(4); %The energy equation
fdot(5)=f(2).*f(4); %The integral energy balance

