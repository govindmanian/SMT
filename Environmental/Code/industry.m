function  [optBeta, RMSEmat, stddev, optdim, rsqall] = industry(X, y, t);

n = size(y,1); %Number of samples
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

for iter=1:t
    iter
    
    %t-fold cross validation. NOT venetian.
    trsize   =   size(X,1)/t;
    cal      =   1:trsize            ;
    val      =   (length(cal) + 1):size(X,1);
    var      =   1:size(X,2)       ;
    Xtr      =   X(cal,var)  ;
    ytr      =   y(cal,:)    ;
    Xte      =   X(val,var)  ;
    yte      =   y(val,:)    ;
    
    Xtrraw = Xtr;
    Xteraw = Xte;
   
    
    %Find the wavelengths in X with the highest correlation to the output
    for cols = 1:size(Xtr,2)
        rsq = corrcoef(ytr,Xtr(:,cols));
        rsq = rsq(2,1)^2;
        rsqall(cols) = rsq;
    end
    
    [~,index] = sort(rsqall, 'descend');
    
    
    %Arbitrary restriction on max dimensionality
    trsize = 10;
    
    if trsize > size(Xtr,1)
        disp('Too many wavelengths - note training block size! Setting max dim to Xtr size')
        trsize = size(Xtr,1)
    end
    
    for dim = 1:trsize
                
        Xtr = Xtrraw(:, index(1:dim));
        Xte = Xteraw(:, index(1:dim));
               
        RMSE = inf;
        
        blr = pinv(Xtr) * ytr;
        yOLS = Xte * blr;
        RMSEraw(dim, iter) = (trace((yte - yOLS)' *(yte - yOLS)) / numel(yte)).^(1/2) ;
        
        if RMSEraw(iter) < RMSE
            RMSE = RMSEraw(dim, iter);
            optBeta = blr;
        end
    end
    
    X=circshift(X,[n2,0]);
    y=circshift(y,n2);
    
end

for row = 1:size(RMSEraw,1)
    RMSEmat(row, 1) = mean(RMSEraw(row,:),2);
    RMSEmat(row, 2) = min(RMSEraw(row,:));
    RMSEmat(row, 3) = max(RMSEraw(row,:));
    stddev(row,1) = std(RMSEraw(row,:),0,2);
end

[~, optdim] = min(RMSEmat(:, 1));