%% Generalised Smooth of Second Order Linear ODE
clear;
clc;

% %Second order linear ODE
x=(0:.05:5)';
n=length(x);

load 'dataerr_0.01.mat'
%load 'dataerr_0.03.mat'
%load 'dataerr_0.05.mat'
%load 'dataerr_0.1.mat'
%plot(x,yerr,'r');

%% Set up the basis with a knot at every point

knots    = x;
norder   = 6;
orderODE = 2;
nbasis   = n + (norder - 2);
rng      = [x(1),x(end)];
basisobj = create_bspline_basis(rng, nbasis, norder, knots);
B        = eval_basis(x, basisobj);

%% Indicate the derivatives to be estimated
% [D0,D1]
Lcoef        = [1,1]';

%% Generalised Smoothoing
argvals=x;  Term=Lcoef; 

for i=1:100
tic;
    [ beta_opt{i}, C_opt{i}, y_hat{i}, var_beta{i}, var_c{i}, pars{i}] = Gen_Pen(argvals, yerr(:,i), basisobj, Term);
time(i)=toc;
end

par = cell2mat(pars);

ind = find(par(1,:)==mode(par(1,:)));

tb = cell2mat(beta_opt);

est_v = cell2mat(var_beta);

%m_SE = cell2mat(SE);

%Optimal beta
display(['Beta ', num2str(mean(tb(:,ind),2)')]);
display(['Var Beta ', num2str(var(tb(:,ind),1,2)')]);
display(['Std Beta ', num2str(std(tb(:,ind),1,2)')]);

%Estimated variance of beta
display(['Est Var Beta ',num2str(mean(est_v(:,ind),2)')]);

%Estimated standard deviation of beta
display(['Est Std Beta ',num2str(sqrt(mean(est_v(:,ind),2)'))]);

%Optimal lambda, SSE, Degress of Freedom
par = cell2mat(pars);
display(['Lambda, SSE, Df ', num2str([mean(par(1,ind)),mean(par(2,ind)),mean(par(3,ind))])]);

%Optimal Sigma Error
display(['Estimated Sigma Error ', num2str(mean(par(2,ind))./(n-mean(par(3,ind))))]);

%Lambda 
tabulate(par(1,:))
%Computation Time
display(['Time ',num2str(mean(time(ind)))]);

%ye = cell2mat(y_hat);

% figure();
% plot(x,yerr(:,ind(1)),'*b');
% hold on;
% plot(x,ye(:,ind(1)),'g');
% plot(x,ye(:,ind(1))+1.96.*mean(m_SE,2),'r');
% plot(x,ye(:,ind(1))-1.96.*mean(m_SE,2),'r');

%save('GenPen0.01.mat')
%save('GenPen0.03.mat')
%save('GenPen0.05.mat')
%save('GenPen0.1.mat')

%save('GenPen0.01_2.mat')
%save('GenPen0.03_2.mat')
%save('GenPen0.05_2.mat')
%save('GenPen0.1_2.mat')



