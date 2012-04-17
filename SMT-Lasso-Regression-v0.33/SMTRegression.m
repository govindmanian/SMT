function bsmt=SMTRegression(y,X,method,option)
%%% with exact cross-validation %%%
 
[n,p]=size(X);

[E,Lambda,SMTArray]=SMTCovarEst(X');
[lambda,I] = sort(diag(Lambda),'descend');
lambda_inv=lambda.^(-1);
E=E(:,I);
Xtilde=X*E*diag(lambda_inv.^(1/2));
btilde=Xtilde'*y/n;

if strcmp(method,'SURE')
%%% SURE %%%
     value=SURE(btilde,n);
     bshrk=Shrk(btilde,value,option);
elseif strcmp(method,'CV')
%%% Cross-validation %%%
t=3;
n1=ceil(n*1/3);
n2=n-(t-1)*n1;
val=0.000:0.01:1;
for iter=1:t
    Xtr=X(1:(t-1)*n1,:);
    ytr=y(1:(t-1)*n1);
    Xte=X((t-1)*n1+1:n,:);
    yte=y((t-1)*n1+1:n);

    Xtetilde=Xte*E*diag(lambda_inv.^(1/2));
    for i=1:size(val,2)
        [bshrk_tmp]=Shrk(btilde,val(i),option);
        snr(t,i)=-10*log10(norm(Xtetilde*bshrk_tmp-yte)^2/n1);
    end
    X=circshift(X,[n2,0]);
    y=circshift(y,n2);
end

avgsnr=mean(snr,1);
[snrmax imax]=max(avgsnr);
value=val(imax);
%figure,plot(val,avgsnr)

bshrk=Shrk(btilde,value,option);


else
   fprintf('Wrong method option\n');
   bshrk=0;
   value=0;
end

bsmt=E*diag(lambda_inv.^(1/2))*bshrk;

