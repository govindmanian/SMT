function blasso=LassoRegression(y,X,method,option)
%%% with exact cross-validation %%%
sigma=1.0; 
[n,p]=size(X);

%%% Cross-validation %%%
t=3;
n1=ceil(n*1/3);
n2=n-(t-1)*n1;
val=5:-1:0;
%val=0:1:5;
for iter=1:t
    Xtr=X(1:(t-1)*n1,:);
    ytr=y(1:(t-1)*n1);
    Xte=X((t-1)*n1+1:n,:);
    yte=y((t-1)*n1+1:n);

    if strcmp(option,'SMTNormalization');
        [E,Lambda,SMTArray]=SMTCovarEst_withoutSORT(Xtr');
    end
    
    for i=1:size(val,2)
   %        [bshrk_tmp]=LassoICD(Xtr,ytr,val(i),10000,option,0, 0);
        if strcmp(option,'SMTNormalization');
            if (i==1 && iter==1)
                [bshrk_tmp]=LassoICD(Xtr,ytr,val(i),10000,option,E, Lambda);
                b_temp=bshrk_tmp;
            else
                [bshrk_tmp]=LassoICD(Xtr,ytr,val(i),10000,option,E, Lambda,b_temp);
                b_temp=bshrk_tmp;
            end
        elseif strcmp(option,'WithNormalization');
            if (i==1 && iter==1)
                [bshrk_tmp]=LassoICD(Xtr,ytr,val(i),10000,option,0, 0);
                b_temp=bshrk_tmp;
            else
                [bshrk_tmp]=LassoICD(Xtr,ytr,val(i),10000,option,0, 0,b_temp);
                b_temp=bshrk_tmp;
            end
        end
        %[bshrk_tmp]=LassoICDTest(Xtr,ytr,val(i),10000,option);
        snr(iter,i)=n*sigma^2/norm(Xte*bshrk_tmp-yte)^2;
    end
    X=circshift(X,[n2,0]);
    y=circshift(y,n2);
end

avgsnr=mean(snr,1);
[snrmax imax]=max(avgsnr);
value=val(imax);
%figure;
%title('For the CV of Lasso');
%plot(avgsnr);

%figure,plot(val,avgsnr)

blasso=LassoICD(X,y,value,10000,option);
%blasso=LassoICDTest(X,y,value,10000,option);
blasso=blasso;

