%%% SMT covariance matrix estimate %%%%%
function [T,Lambda,SMTArray]=SMTCovarEst(Y,K)

[N,M]=size(Y);

if nargin == 1
   index=randperm(M);
   Y=Y(:,index);
   M1=ceil(M*1/3);

   
   %NumRotations = ceil(3*log2(M)*N);
   
   NumRotations = ceil(log2(M)*N);

   Rhat = Y(:,1:2*M1)*Y(:,1:2*M1)'/(2*M1);
   [T,lambda,SMTArray,likelihood1] = ComputeSMT(Rhat,NumRotations,Y(:,2*M1+1:M));

   Y=circshift(Y,[0,M-2*M1]);
   Rhat = Y(:,1:2*M1)*Y(:,1:2*M1)'/(2*M1);
   [T,lambda,SMTArray,likelihood2] = ComputeSMT(Rhat,NumRotations,Y(:,2*M1+1:M));

   Y=circshift(Y,[0,M-2*M1]);
   Rhat = Y(:,1:2*M1)*Y(:,1:2*M1)'/(2*M1);
   [T,lambda,SMTArray,likelihood3] = ComputeSMT(Rhat,NumRotations,Y(:,2*M1+1:M));

   likelihood=likelihood1+likelihood2+likelihood3;
   [valmax,K]=max(likelihood);
   %figure,set(gca,'fontsize',18)
   %plot(likelihood,'LineWidth',2)
   %xlabel('# of SMT rotations')
   %ylabel('Average log-likelihood')
   if K==NumRotations
       disp('It may require more rotations!');
   end
end

K;

Rhat = (1/M)*Y*Y';
[T,lambda,SMTArray,likelihood] = ComputeSMT(Rhat,K);
[lambda,I] = sort(lambda);
T=T(:,I);
Lambda=diag(lambda);

