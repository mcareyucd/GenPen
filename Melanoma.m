%% Generalised Smooth of Melanoma data
clear;
clc;

%% Load the Melanoma data
tempmat = load('melanoma.dat');

x  = tempmat(:,2);
y  = tempmat(:,3);

n=length(y);

plot(x,y)

%% Set up the basis with a knot at every point

knots    = x;
norder   = 6;
orderODE = 4;
nbasis   = n + (norder - 2);
rng      = [x(1),x(end)];
basisobj = create_bspline_basis(rng, nbasis, norder, knots);
B        = eval_basis(x, basisobj);

%% Indicate the derivatives to be estimated
% [D0,D1,D2,D3]
Lcoef        = [0,0,1,0]';

%% Generalised Smoothoing
 argvals=x;  Term=Lcoef;
tic;
[ beta_opt, C_opt, y_hat, var_beta, var_c, pars] = Gen_Pen(argvals, y, basisobj, Term);
toc;
%[ beta_opt_2,C_opt_2,lambda_opt_2,y_hat_2,method_1,method_2, SumSE_2, df_2 ] = Gen_Pen_2(argvals, y, basisobj, Term);
%Optimal beta
display(beta_opt);
%Variance of beta
display((var_beta));
display(sqrt(var_beta));

%Optimal lambda
display(pars);

%Plot
figure();
plot(x,y,'b.')
hold on;
plot(x,y_hat,'k-')


% %%
% % Heckman and Ramsay (2000) method and Cao and Ramsay (2007) method
% 
% %  Last modified  26 July 2006 by Jim Ramsay: code http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/
% 
% %  -----------------------------------------------------------------------
% %                      Melanoma Incidence data
% %  -----------------------------------------------------------------------
% 
% %  input the data  
% 
% tempmat = load('melanoma.dat');
% 
% year  = tempmat(:,2);
% mela  = tempmat(:,3);
% nyear = length(year);
% rng   = [min(year),max(year)];
% 
% %  -----------------------------------------------------------------------
% %             smooth data using B-splines
% %  -----------------------------------------------------------------------
% 
% %  set up the basis with a knot at every year
% 
% knots    = year';
% nbasis   = nyear + 4;
% norder   = 6;
% basisobj = create_bspline_basis(rng, nbasis, norder, knots);
% 
% %  set up operator to remove sinusoid plus trend
% 
% %Heckman and Ramsay
% %omega  = .3365;
% %Cao and Ramsay
% omega  = .4106;
% 
% Lbasis = create_constant_basis(rng);
% Lcoef  = [0, 0, omega, 0];
% wfd          = fd(Lcoef, Lbasis);
% wfdcell      = fd2cell(wfd);     % convert the FD object to a cell object
% melaLfd      = Lfd(4, wfdcell);  %  define the operator object
% %Heckman and Ramsay
% %lambda       = .36.*10.^5;
% %Cao and Ramsay
% lambda       = 2.4.*10.^4;
% melafdPar    = fdPar(basisobj, melaLfd, lambda);
% 
% %  smooth the data
% 
% melafd = smooth_basis(year, mela, melafdPar);
% 
% %yhat_HR = eval_fd(year,melafd,0);
% yhat_RamCao = eval_fd(year,melafd,0);