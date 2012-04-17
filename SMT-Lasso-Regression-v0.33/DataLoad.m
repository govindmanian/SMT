function [X, y, Xte, yte, Rw, t]=DataLoad(filename, n_training_sample, option);
disp('%%%%%%%%%%%%%%%%%%%%%%DATA LOADING&GENERATING%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('                 ');
%The training data and the testing data are generated seprately by
%X=y*t+W', where 
%size(t)=px1, t is a deterministic but unknown signal;
%size(y)=nx1, y is a vector of n independent Gaussian random variables to
%be estimated;
%size(W)=nxp, each row of W is an independent p-dimensional Gaussian random
%vector of correlated noise or clutter. 
%Notice: Rw, which is the covariance matrix of W, is computed from the real
%byperspectral data
%More descprition about the experiment data can be found in the
%corresponding paper:
%Guangzhi Cao, Yandong Guo, and Charles A. Bouman, 
%``High Dimensional Regression using the Sparse Matrix Transform (SMT),'' 
%in the Proceedings of the {\em International Conference on Acoustic, Speech, and Signal Processing (ICASSP)}, 
%March 14-19, 2010.

%Configuration
load (filename)
n=n_training_sample;  % 50, 100, or 200
s=3;

W=pixels;

[p N]=size(W);

%% compute 'true' R
avgW=mean(W,2);
W=W-avgW*ones(1,N);
Rw=1/N*W*W';   
[Ew,Lambdaw] = eig(Rw);
[lambdaw,I] = sort(diag(Lambdaw),'descend');
Ew=Ew(:,I);
Lambdaw=diag(lambdaw);

if strcmp(option,'random');
    t=randn(1,p); %t can be random or one of the column of Ew
elseif strcmp(option,'eigenv');
    t=Ew(:,170)';
else
    fprintf('Wrong method option\n');
    return;
end

t=s*t/norm(t);

%training data
sigma=1;
W = randn([n,p])*sqrt(Lambdaw)*Ew';   %% Generate Gaussian samples from covariance R
y=randn(n,1)*sigma;
X=y*t+W;

%testing data
nt=300;
yte=randn(nt,1);
Wte=randn([nt,p])*sqrt(Lambdaw)*Ew';
Xte=yte*t+Wte;