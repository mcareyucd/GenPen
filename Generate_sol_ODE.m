clc;
clear;

% Load Simulated data
x=(0:.05:5)';
n = length(x);

%coefficient in front of highest order term in ODE (D2)
a=1;
%coefficient in front of second highest order term in ODE (D1)
b=2;
%coefficient in front of lowest order term in ODE (D0)
c=10;

%discriminant
dis=b.^2-4*a*c;

%roots of auxiliary equation
r1=(-b+sqrt(dis))/(2*a);
r2=(-b-sqrt(dis))/(2*a);

randn('seed',0);

if(dis>0)
    %the general solution
    y=exp(r1.*x)+exp(r2.*x);
end
if(dis==0)
    %the general solution
    y=exp(r1.*x)+x.*exp(r1.*x);
end
if(dis<0)
    %the general solution
    alpha=-b/(2*a);
    beta1=sqrt(-dis)/(2*a);
    y=exp(alpha.*x).*(cos(beta1.*x)+sin(beta1.*x));
    dy=exp(alpha.*x).*(alpha.*cos(beta1.*x)+alpha.*sin(beta1.*x)-sin(beta1.*x).*beta1+cos(beta1.*x).*beta1);
end

rng('default');

err = randn(n,1000);
yerr=repmat(y,[1,1000])+0.01.*err;
%dyerr=repmat(dy,[1,1000])+0.01.*err;

save('dataerr_0.01.mat', 'yerr');

err = randn(n,1000);
yerr=repmat(y,[1,1000])+0.02.*err;
%dyerr=repmat(dy,[1,1000])+0.02.*err;

save('dataerr_0.02.mat', 'yerr');

err = randn(n,1000);
yerr=repmat(y,[1,1000])+0.03.*err;
%dyerr=repmat(dy,[1,1000])+0.03.*err;

save('dataerr_0.03.mat', 'yerr');

err = randn(n,1000);
yerr=repmat(y,[1,1000])+0.05.*err;
%dyerr=repmat(dy,[1,1000])+0.05.*err;

save('dataerr_0.05.mat', 'yerr');

err = randn(n,1000);
yerr=repmat(y,[1,1000])+0.1.*err;
%dyerr=repmat(dy,[1,1000])+0.1.*err;

save('dataerr_0.1.mat', 'yerr');