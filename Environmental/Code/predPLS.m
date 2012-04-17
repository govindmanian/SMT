%%Predict y based on PLS model

function Ypred = predPLS(model,X,k) ;

W       =   model.W(:,1:k)     ;
Q       =   model.Q(:,1:k)     ;
beta    =   W*Q'                ;
Ypred   =   X*beta  ;