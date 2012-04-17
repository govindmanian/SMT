%%Ridge Regression

function [optBeta, RMSEmat, minLambda, edof, minedof, R2,ridgestd] = Ridge(X, y, t)

[n, m] = size(X); %Number of samples

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
    
    mtr = size(Xtr,2);
    mte = size(Xte,2);
    mbs(1,iter) = mte;

    count = 1;
    
    %Set regularization
    meta_min    =   -1   ;
    meta_max    =   7   ;
    meta_res    =   30   ;
    domain = 10.^(meta_min:abs(meta_max/meta_res):meta_max)    ;


    for lambda = domain
       
        kernel = (Xtr' * Xtr + lambda * eye(mtr,mtr)) \ Xtr'; 
        beta = kernel * ytr; %Calculate beta from training set    
                
        yhat = Xte * beta;

        RMSE = (trace((yte - yhat)' *(yte - yhat)) / numel(yte)).^(1/2);

        RMSEraw(count, iter) = RMSE;
        
        edofraw(count, iter) = trace(Xtr * kernel);
                
        count = count + 1;
        
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

%Calculate R2 for all blocks
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
    ridgestd(row,1) = std(RMSEraw(row,:),0,2);
    edof(row, 1) = mean(edofraw(row,:),2);
    edof(row, 2) = min(edofraw(row,:));
    edof(row, 3) = max(edofraw(row,:));
end



%Find the model parameters by fitting over THE ENTIRE DATA SET
[~, index] = min(RMSEmat(:, 1));
minedof = edof(index,1);
optBeta = (X' * X + domain(index) * eye(m,m)) \ X' * y; %Calculate beta from training set

minLambda = domain(index);
end