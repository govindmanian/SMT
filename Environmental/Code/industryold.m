function  [optBeta, RMSEmat, stddev, optdim, rsqall] = industry(X, y, t);

n = size(y,1); %Number of samples
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

for iter=1:t
    iter
    
    %t-fold cross validation
    trsize   =   size(X,1)/t;
    cal      =   1:trsize            ;
    val      =   (length(cal) + 1):size(X,1);
    var      =   1:size(X,2)       ;
    Xtr      =   X(cal,var)  ;
    ytr      =   y(cal,:)    ;
    Xte      =   X(val,var)  ;
    yte      =   y(val,:)    ;
    
    %Find the wavelengths in X with the highest correlation to the output
    for cols = 1:size(Xtr,2)
        rsq = corrcoef(ytr,Xtr(:,cols));
        rsq = rsq(2,1)^2;
        rsqall(cols) = rsq;
    end
    
    %Find the wavelength where the max is at
    [~, index] = max(rsqall);
    
    %Arbitrary restriction on dimensionality
    trsize = 10;
    
    if index - trsize/2 <= 1
        Xtr = Xtr(:, 1:trsize);
        Xte = Xte(:, 1:trsize);
    elseif index + trsize/2 > size(Xtr,2)
        Xtr = Xtr(:, end - trsize:end);
        Xte = Xte(:, end - trsize:end);        
    else
        Xtr = Xtr(:, index - trsize/2:index + trsize/2);
        Xte = Xte(:, index - trsize/2:index + trsize/2);
    end 
    
    RMSE = inf;
    
    blr = pinv(Xtr) * ytr;
    yOLS = Xte * blr;
    RMSEraw(1, iter) = (trace((yte - yOLS)' *(yte - yOLS)) / numel(yte)).^(1/2) ;
    
    if RMSEraw(iter) < RMSE
        RMSE = RMSEraw(iter);
        optBeta = blr;
        optdim = index;
    end
    
    X=circshift(X,[n2,0]);
    y=circshift(y,n2);
    
end

RMSEmat = RMSEraw;
stddev = std(RMSEmat);