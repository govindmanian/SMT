%The software can estimate the linear regression coefficients by 
%traditional lasso, SMT-Based lasso, SMT-shrinkage, etc.
%The descriptions of these algorithms can be found in the paper
%Guangzhi Cao, Yandong Guo, and Charles A. Bouman, 
%``High Dimensional Regression using the Sparse Matrix Transform (SMT),'' 
%in the Proceedings of the {\em International Conference on Acoustic, Speech, and Signal Processing (ICASSP)}, 
%March 14-19, 2010.

clear all
close all

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('STM-Lasso-Regression-v0.31');
disp('                 ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('                 ');

%Data loading process:
%Load the data used in the estimation.
%X, y are the training pairs, which are used to estimate the regression
%vector;
%Xte and yte are the testing pairs, which are used to evaluate the
%performance of the regression vector got by the training data.
%Rw and t are used to calculate the covariance matrix 'Rx' for X, which are
%just for testing;
%Users can apply our algorithms to their own training and testing data by
%subsituting this following function. 
filename='grass'; % can be 'grass' or 'water';
option='eigenv';  % can be 'random' or 'eigenv';
[X, y, Xte, yte, Rw, t]=DataLoad(filename, 100, option);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Training process%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%TRAINING PROCESS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('                 ');

[n,p]=size(X);
[nt, p]=size(Xte);
sigma=1;

%%%%%%%%%%%%%%%%%%%%%%  zero estimator  %%%%%%%%%%%%%%%%%%%
disp('Set all the coefficients to be zero');
b0=0;
disp('Estimation accomplished!');
disp('                 ');

%%%%%%%%%%%%%%%%%%%%%%  Traditional Lasso  %%%%%%%%%%%%%%%%%%%
%This method is time consuming, so it is turned off by default
RunLasso = 0;
if RunLasso
  tic;
  LargerDim='p';
  method='CV'; % 'CV' or 'SURE'
  option='WithNormalization';
  disp('Traditional Lasso Regression with Cross Validation');
  disp('Commencing...');
  blasso=LassoRegression(y,X,method,option);
  toc;
  disp('Estimation accomplished!');
  disp('                 ');
end

%%%%%%%%%%%%%%%%%%%%%%  SMT-Lasso  %%%%%%%%%%%%%%%%%%%
tic;
method='CV'; % 'CV' or 'SURE'
option='SMTNormalization';
disp('SMT Based Lasso Regression with cross validation');
disp('Commencing...');
blasso_SMT=LassoRegression(y,X,method,option);
toc;
disp('Estimation accomplished!');
disp('                 ');

%%%%%%%%%%%%%%%%%%%%%%  SMT-shrinkage  %%%%%%%%%%%%%%%%%%%
method='CV'; % 'CV' or 'SURE'
option='soft';
disp('SMT shrinkage Regression');
disp('Commencing...');
bsmt=SMTRegression(y,X,method,option);
disp('Estimation accomplished!');
disp('                 ');

%%%%%%%%%%%%%%%%%%%%%%  SMT-subset  %%%%%%%%%%%%%%%%%%%
method='CV'; % 'CV' or 'SURE'
option='hard';
disp('SMT subset Regression');
disp('Commencing...');
bsmt_S=SMTRegression(y,X,method,option);
disp('Estimation accomplished!');
disp('                 ');

%%%%%%%%%%%%%%%%%%%%%%  SMT-shrinkage with SURE  %%%%%%%%%%%%%%%%%%%
method='SURE'; % 'CV' or 'SURE'
option='soft';
disp('SMT Shrinkage regression with SURE');
disp('Commencing...');
bsmt_SURE=SMTRegression(y,X,method,option);
disp('Estimation accomplished!');
disp('                 ');


%%%%%%%%%%%%%%%%%%%%  OLS without Rx %%%%%%%%%%%%%%%%%%%
disp('Original least square estimation without known the true covariance');
disp('Commencing...');
LargerDim='n';
method='CV'; % 'CV' or 'SURE'
%if n>p;
    bols=LinearRegression(y,X,LargerDim,'N');
%else p>n;
%    bols=zeros(p,1);
%end
disp('Estimation accomplished!');
disp('                 ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Testing process%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%TESTING PROCESS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('                 ');

snr(1)=nt*sigma^2/norm(b0-yte)^2;
if RunLasso
 snr(2)=nt*sigma^2/norm(Xte*blasso-yte)^2;
end
snr(3)=nt*sigma^2/norm(Xte*blasso_SMT-yte)^2;
snr(4)=nt*sigma^2/norm(Xte*bsmt-yte)^2;
snr(5)=nt*sigma^2/norm(Xte*bsmt_S-yte)^2;
snr(6)=nt*sigma^2/norm(Xte*bsmt_SURE-yte)^2;
snr(7)=nt*sigma^2/norm(Xte*bols-yte)^2;

disp('SNR for different estimation methods');
%SNR is defined by one over mean square error;
fprintf('%20s %20s %20s %20s\n','Zeros vector','Traditional Lasso','SMT-Lasso','SMT-shrinkage');
if RunLasso
    fprintf('%20.4f %20.4f %20.4f %20.4f\n',snr(1),snr(2),snr(3),snr(4));
end
fprintf('%20.4f %20s %20.4f %20.4f\n',snr(1),'NAL',snr(3),snr(4));
fprintf('%20s %20s %20s \n','SMT-subset','SMT-SURE','Least Square');
fprintf('%20.4f %20.4f %20.4f\n',snr(5),snr(6),snr(7));

save ./results/snr-grass-100.mat snr
