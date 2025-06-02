function [lambda,eigenvec]=slsolve(varargin)
%This function solves the Sturm-Liouville eigenvalue problem given by:
%
%  [p(x) y']' - q(x) y + lambda w(x) y = 0
%
% subject to the boundary conditions:
%
%  bc(1) y(0) + bc(2) y'(0) = 0 
%
%  bc(3) y(1) + bc(4) y'(1) = 0
%
% over the domain 0< x <1.
%
% The function is called by the command:
%
% [lambda,eigenvec]=slsolve('pfun','qfun','wfun',bc,n);
%
% The function call requires that you provide the function names (or
% handles if you are using the anonymous function utility) for the
% functions p, q, and w.  These functions must be able to handle an array
% of values.  You also provide the boundary coefficients in the array bc.
% In addition, you may specify the degree of discretization n.  Its default
% value is 50.  The matrices which are generated are of size (n+1,n+1).
% The function returns the eigenvalues in the array lambda (sorted by size
% in ascending order) and the matrix eigenvec which contains the
% corresponding eigenfunctions.  The eigenfunctions are all normalized by
% their maximum value over the domain 0<x<1.
%
% A last note on error:  The code uses second order derivative
% approximations, so the error in the eigenvalues and eigenvectors will be
% of O(1/n^2).  In general, the first few eigenvalues will be reliable, but
% the accuracy will deteriorate as you look at the higher eigenvalues, with
% the last few being meaningless.

p=varargin{1};
q=varargin{2};
w=varargin{3};
bc=varargin{4};

if nargin<5;n=50;else;n=varargin{5};end

h=1/n; %set discretization
x=[0:h:1]'; %this is the array of x values

%Now we set up the arrays used in making the matrix A:
pp=zeros(1,n+1);
pm=zeros(1,n+1);
ww=zeros(1,n+1);
qq=zeros(1,n+1);
for i=2:n
 pp(i)=feval(p,x(i)+h/2);
 pm(i)=feval(p,x(i)-h/2);
 ww(i)=feval(w,x(i));
 qq(i)=feval(q,x(i));
end

%The matrix W is easy:
weight=-diag(ww);

%The matrix A is a bit more complex.  First we do
%the main diagonal:
a=diag(-pp-pm-qq*h^2);
%and then the super and sub diagonals:
a=a+diag(pp(1:n),1);
a=a+diag(pm(2:n+1),-1);

%Finally, we divide by h^2:
a=a/h^2;

%And now for the boundary conditions.  First at the left edge:
a(1,1)=bc(1)-bc(2)*1.5/h;
a(1,2)=bc(2)*2/h;
a(1,3)=-bc(2)/2/h;

%and at the right edge:
a(n+1,n+1)=bc(3)+bc(4)*1.5/h;
a(n+1,n)=-bc(4)*2/h;
a(n+1,n-1)=bc(4)/2/h;

%Now we are ready to calculate the eigenvalues:
[v,d]=eig(a,weight);

%The number of eigenvalues and vectors will be less
%than the size of A and W, thus:
evals=diag(d);
i=find(isfinite(evals)); %The matlab 7 form of finite!
evals=evals(i);
evecs=v(:,i);

%Now we sort the eigenvalues and eigenvectors according
%to the size of the eigenvalues:
[~,i]=sort(abs(real(evals)));
lambda=evals(i);
evecs=evecs(:,i);

%and finally, we normalize the eigenvectors by their
%maximum value.
eigenvec=zeros(size(evecs));
for j=1:length(lambda)
 eigenvec(:,j)=evecs(:,j)/norm(evecs(:,j),inf)/sign(evecs(2,j));
end