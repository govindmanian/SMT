function [optBeta, RMSEmat, R2] = OLS(X, y, t)


n = size(y,1); %Number of samples
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

RMSE = inf;


for iter=1:t
    iter

    cte = 1;
    ctr = 1;
    %Use venetian blinds for testing and training sets
    for ele = 1:size(y,1)

        if mod(ele,t) == 1 %If element is mod t, put in training set (so tr:te is 1:(t-1))
            ytr(ctr,:) = y(ele,:);
            Xtr(ctr,:) = X(ele,:);
            ctr = ctr + 1;

        else
            yte(cte,:) = y(ele,:);
            Xte(cte,:) = X(ele,:);
            cte = cte + 1;

        end
    end
    
    mte = size(Xte,2);
    mbs(1,iter) = mte;

    
    blr = pinv(Xtr) * ytr;
    yOLS = Xte * blr;
    RMSEraw(iter) = (trace((yte - yOLS)' *(yte - yOLS)) / numel(yte)).^(1/2) ;
    
    if RMSEraw(iter) < RMSE
        RMSE = RMSEraw(iter);
        optBeta = blr;
    end
    
    vary(1,iter) = var(yte);
    
    %Using venetian blinds, so just shift the elements over by 1
    X=circshift(X,[1,0]);
    y=circshift(y,1);
end

RMSEmat = RMSEraw;

%Calculate R2 for all
denom = sum(mbs .* vary); %Denominator

for column = 1:size(mbs,2) %Duplicate testing block size to do elementwise multiplication
    mbs(1:size(RMSEraw,1), column) = mbs(1,column);
end

RMSEave2 = mbs .* RMSEraw.^2; 

for row = 1:size(RMSEave2,1) %Average over iterations
    numer(row,1) = mean(RMSEave2(row,:),2);
end

R2 = 1 - numer/denom;