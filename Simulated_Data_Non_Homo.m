% %% Generalised Smoothoing for simation of non-homogenous ODE
clc;
clear;

%Set up the B-spline basis function
x = (0:.5:4.5)';
n=length(x);
knots    = x;
norder   = 4;
nbasis   = n + (norder - 2);

rng      = [x(1),x(end)];
basisobj = create_bspline_basis(rng, nbasis, norder, knots);
Term = [1,1]';

%Solution of the non-homogenous ODE
n=length(x);
y = zeros(n,1);
for i = 1:length(x)
y(i)=(47300/357)*exp(-(133/200)*x(i))*sin((3/200)*sqrt(119)*x(i))*sqrt(119)+100*exp(-(133/200)*x(i))*cos((3/200)*sqrt(119)*x(i))+...
(50/51)*sqrt(119)*heaviside(x(i)-1/10)*exp(133/2000-(133/200)*x(i))*sin((1/2000*(-3+30*x(i)))*sqrt(119));
end

plot(x,y)

%add error to generate 1000 samples
err = randn(n,1000);
yerr=repmat(y,[1,1000])+0.01.*err;

beta_opt  = cell(1000,1);
C_opt     = cell(1000,1);
y_hat     = cell(1000,1);
var_beta  = cell(1000,1);
pars      = cell(1000,1);

time = zeros(100,1);

for i=1:1000
tic;
    [beta_opt{i}, C_opt{i}, y_hat{i}, var_beta{i}, pars{i}] = Gen_Pen_Non(x, yerr(:,i), basisobj, Term);
time(i)=toc;
end

par = reshape(cell2mat(pars),3,1000);

ind = find(par(1,:)==mode(par(1,:)));

tb = reshape(cell2mat(beta_opt),2,1000);

est_v = reshape(cell2mat(var_beta),2,1000);

%Optimal beta
display(['Beta ',     num2str(mean(tb(:,ind),2)')]);
display(['Var Beta ', num2str(var(tb(:,ind),1,2)')]);
display(['Std Beta ', num2str(std(tb(:,ind),1,2)')]);

% %Estimated variance of beta
display(['Est Var Beta ',num2str(mean(est_v(:,ind),2)')]);
 
% %Estimated standard deviation of beta
display(['Est Std Beta ',num2str(sqrt(mean(est_v(:,ind),2)'))]);
  
% %Optimal lambda, SSE, Degress of Freedom
display(['Lambda, SSE, Df ', num2str([mean(par(1,ind)),mean(par(2,ind)),mean(par(3,ind))])]);
% 
% %Optimal Sigma Error
display(['Estimated Sigma Error ', num2str(mean(par(2,ind))./(n-mean(par(3,ind))))]);
 
%Time
display(['Time',num2str(mean(time))]);

