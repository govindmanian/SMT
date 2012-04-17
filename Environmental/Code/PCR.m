%SVD decomposition PCR algorithm

function [optBeta, RMSEmat, optdim, R2, pcrstd] = PCR(X, y, t)

n = size(y,1);
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

minRMSE = inf;

for iter=1:t
    iter
    
    cte = 1;
    ctr = 1;
    
    %Use venetian blinds for cross validation
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
    
    
    [U,S,P] = svd(Xtr,'econ');
    lambda = diag(S); %/ntr;   % eigenvalues
    T = U*diag(lambda);
    
    
    pc = 15 % select 5 PC's (model dimension)
    
    if pc > size(Xtr,1)
        disp('Too many pcs - note training block size! Setting max dim to Xtr size')
        pc = size(Xtr,1)
    end
    
    %     Preallocate/initialize
    
    for dim = 1:pc
        dim
        Ttr = T(:,1:dim);
        Ptr = P(:,1:dim);
        
        beta = (Ttr' * Ttr) \ Ttr' * ytr;
        betaX = Ptr * beta;
        yhat = Xte * betaX;
        
        RMSE = (trace((yte - yhat)' *(yte - yhat)) / numel(yte)).^(1/2) ;
        RMSEraw(dim, iter) = RMSE;
        
        if RMSE < minRMSE
            minRMSE = RMSE;
            kyhat = yhat;
            kyte = yte;
        end
        
        vary(1,iter) = var(yte);
        
    end
    
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

%Generate mean, max, min for each component
for row = 1:size(RMSEraw,1)
    RMSEmat(row, 1) = mean(RMSEraw(row,:),2);
    RMSEmat(row, 2) = min(RMSEraw(row,:));
    RMSEmat(row, 3) = max(RMSEraw(row,:));
    pcrstd(row,1) = std(RMSEraw(row,:),0,2);
end

%Find the model parameters by fitting over THE ENTIRE DATA SET
[U,S,P] = svd(X,'econ');
lambda = diag(S);
T = U*diag(lambda);
[~, optdim] = min(RMSEmat(:, 1));
optBeta = P(:,1:optdim) * ((T(:,1:optdim)' * T(:,1:optdim)) \ T(:,1:optdim)' * y);

end