function blr=LinearRegression(y,X,LargerDim,SMT)
%%% with exact cross-validation %%%
 
[n,p]=size(X);

if strcmp(SMT,'N')
    if strcmp(LargerDim,'p')
        blr=X'*(inv(X*X'))*y;
    elseif strcmp(LargerDim,'n')
        blr=(inv(X'*X))*X'*y;
    end
elseif strcmp(SMT,'Y')
    [E,Lambda,SMTArray]=SMTCovarEst(X');
    [lambda,I] = sort(diag(Lambda),'descend');
    lambda_inv=lambda.^(-1);
    E=E(:,I);
    Xtilde=X*E*diag(lambda_inv.^(1/2));
    if strcmp(LargerDim,'p')
        btilde=Xtilde'*(inv(Xtilde*Xtilde'))*y;
    elseif strcmp(LargerDim,'n')
        btilde=(inv(Xtilde'*Xtilde))*Xtilde'*y;
    end
    blr=E*diag(lambda_inv.^(1/2))*btilde;
end
