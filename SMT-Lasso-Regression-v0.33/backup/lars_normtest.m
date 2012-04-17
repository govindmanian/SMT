function b_lars=lars_normtest(X,y,option);
[n,P]=size(X);
   eps=1e-5;
   if strcmp(option,'WithNormalization');
       mu=mean(X);
       sigma2 = std(X);
       ndx = find(sigma2 < eps);
       sigma2(ndx) = 1;
       E=zeros(P,P);
       for i=1:P;
            E(i,i)=1/sigma2(i);
       end
       X = X *E;
       X(:,P)=ones(n,1);
   elseif strcmp(option,'SMTNormalization');
        %disp('SMT based');
        [E,Lambda,SMTArray]=SMTCovarEst_withoutSORT(X');
            SMTlambda=diag(Lambda);
            SMTlambda_inv=SMTlambda.^(-1);
            Xtilde=X;
            Xtilde=Xtilde*E*diag(SMTlambda_inv.^(1/2));
            X=Xtilde;
            %b=E*diag(SMTlambda_inv.^(1/2))*b; %It seems that to use this multiply directly is even faster...
   else
       disp('Warning: the option has to be either "WithNormalization" or "SMTNormalization"');
   end
    %Initialize the regression vector beta;
[s_opt, b_opt, res_mean, res_std] = crossvalidate(@lars, 10, 1000, X, y, 'lars', 0, 0, [], 0);
    b=b_opt;
    b=b';
    
    if strcmp(option,'WithNormalization');
        b_lars=E*b;
    elseif strcmp(option,'SMTNormalization');
        b_lars=E*diag(SMTlambda_inv.^(1/2))*b; %It seems that to use this multiply directly is even faster...
    else
       disp('Warning: the option has to be either "WithNormalization" or "SMTNormalization"');
    end





