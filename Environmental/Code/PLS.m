%PLS1 algorithm for PLS

function [optBeta, RMSEmat, index, R2, plsstd] = PLS(X, y, t)

n = size(y,1);
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

minRMSE = inf;

for iter=1:t
    iter


    cte = 1;
    ctr = 1;
    %Use venetian blinds
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

    
    kmax = 10; % establish PLS models until this max number of components (latent variables)
    model   =   makePLS(Xtr,ytr,kmax)   ;
    
    for k = 1:kmax   ; % select dimension       
        yhat = predPLS(model, Xte, k);
        RMSE = (trace((yte - yhat)' *(yte - yhat)) / numel(yte)).^(1/2) ;
        
        RMSEraw(k, iter) = RMSE;
        
        if RMSE < minRMSE
            minRMSE = RMSE;
            kyhat = yhat;
            kyte = yte;
        end
    end
    
    vary(1,iter) = var(yte);


    %Using venetian blinds, so just shift the elements over by 1
    X=circshift(X,[1,0]);
    y=circshift(y,1);
end

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


for row = 1:size(RMSEraw,1)
    RMSEmat(row, 1) = mean(RMSEraw(row,:),2);
    RMSEmat(row, 2) = min(RMSEraw(row,:));
    RMSEmat(row, 3) = max(RMSEraw(row,:));
    plsstd(row,1) = std(RMSEraw(row,:),0,2);
end

%Find the model parameters by fitting over THE ENTIRE DATA SET
[~, index] = min(RMSEmat(:, 1));
model   =   makePLS(X,y,index)   ;
optBeta = model.W(:,1:index) * model.Q(:,1:index)';

end