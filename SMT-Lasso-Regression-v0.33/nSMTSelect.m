%function   [T,lambda,SMTArray,likelihood] = ComputeSMT(R,NumRotations,Y)
% R - Input convariance matrix
% NumRotations - number of rotations
% Y - the sample data
%
% SMTArray contains the sequence of SMT rotions. So that
%           SMTArray(i,1) = first index
%           SMTArray(i,2) = second index
%           SMTArray(i,3) = rotation angle
%
% T is the product of the SMT rotations, so that
%           D = T*R*T^t;
% where D is approximately diagonal
% and       lambda = diag(D)

%%%This is a version of SMT that does dimension reduction based on the size
%%%of the eigenvalues

clc;
clear;
close all;

dataset = 1;

if dataset == 1
    load('spectra.mat');
    Xraw = DataSet.Spectra;
    Yraw = DataSet.Conc;
    
    %Find sizes
    n = size(Xraw,1); %Number of samples
    X = Xraw - ones(n,1) * mean(Xraw);
    y = Yraw - ones(n,1) * mean(Yraw);
elseif dataset == 2
    load('INORFULL.mat');
    Xraw = DATA;
    Yraw = CONC;
    n = size(Xraw,1); %Number of samples
    X = Xraw - ones(n,1) * mean(Xraw);
    y = Yraw - ones(n,1) * mean(Yraw)./(ones(n,1) * std(Yraw,0,1));
end

t = 1;
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

initialize = 1;

for iter=1:t
    iter
    
    trsize   =   100;%26;
    cal      =   1:trsize       ;
    val      =   (length(cal) + 1):size(X,1);
    var      =   1:size(X,2)       ;
    Xtr      =   X(cal,var)  ;
    ytr      =   y(cal,:)    ;
    Xte      =   X(val,var)  ;
    yte      =   y(val,:)    ;
    
    [ntr, mtr] = size(Xtr);
    [nte, mte] = size(Xte);
    
    
    NumRotations = ceil(log2(mtr)*ntr);
    
    R=1/mtr*Xtr'*Xtr;
    
    T = eye(size(R));
    N = max(size(R));
    F = Criteria(R);
    [f,Irow]=max(F);
    
    
    for iterations=1:NumRotations
        iterations
        
        %%%% fast search %%%%%%%
        [val,Icol]=max(f);
        i=Irow(Icol);
        j=Icol;
        
        if(i>j)
            ii = i; i=j; j=ii;
        end
        
        rii = R(i,i);
        rij = R(i,j);
        rjj = R(j,j);
        
        [a,b,theta] = rotation(rii,rjj,rij);
        
        %%% fast implementation %%%%
        % R=A'*R*A
        ci = R(:,i)*cos(theta)+R(:,j)*sin(theta);   % R*A
        cj = -R(:,i)*sin(theta)+R(:,j)*cos(theta);
        R(:,i) = ci;
        R(:,j) = cj;
        
        ri = R(i,:)*cos(theta) + R(j,:)*sin(theta); % A'*R*A
        rj = -R(i,:)*sin(theta) + R(j,:)*cos(theta);
        R(i,:) = ri;
        R(j,:) = rj;
        
        
        % T=T*A
        ri = T(:,i)*cos(theta)+T(:,j)*sin(theta);   % T*A
        rj = -T(:,i)*sin(theta)+T(:,j)*cos(theta);
        T(:,i) = ri;
        T(:,j) = rj;
        
        
        %%%%%% update the F matrix %%%%%%
        lambda=diag(R);
        F(:,i)=R(:,i).^2 ./ (R(i,i)*lambda);
        F(:,j)=R(:,j).^2 ./ (R(j,j)*lambda);
        F(i,:)=F(:,i)';
        F(j,:)=F(:,j)';
        F(i,i)=0;
        F(j,j)=0;
        
        %%%%%% update the f and Irow vectors %%%%%%%%
        for z=1:N
            if Irow(z)==i || Irow(z)==j || z==i || z==j
                [f(z),Irow(z)]=max(F(:,z));
            else
                if f(z)<F(i,z)
                    f(z)=F(i,z);
                    Irow(z)=i;
                end
                if f(z)<F(j,z)
                    f(z)=F(j,z);
                    Irow(z)=j;
                end
            end
        end
        
        lambda=diag(R);
        
        [Lambda,I] = sort(lambda);
        E = T(:,I);
        Lambda = diag(Lambda);
        
        Lambda = diag(Lambda);
        [Lambda,I] = sort(Lambda,'descend');
        E = E(:,I);
        
        itercount = 1;
        if mod(iterations,NumRotations) == 0
            
            cutcount = 1;
            
            for index = 1:10
                cutoff              =   Lambda(index);
                sel                 =   find(Lambda>cutoff,1,'last');
                selE                =   E(:,1:sel)                  ;
                selLambda           =   Lambda(1:sel)               ;
                
                Lambda_inv = selLambda.^(-1);
                P = selE*diag(Lambda_inv.^(1/2));
                Ttr = Xtr*P;
                Tte = Xte*P;
                btilde = pinv(Ttr) * ytr;
                
                BetaSMT = P*btilde;
                
                yhat = Tte*btilde;
                RMSE = (trace((yte - yhat)' *(yte - yhat)) / numel(yte)).^(1/2);
                RMSEmat(itercount, index, iter) = RMSE;
                
                if initialize == 1
                    minRMSE = RMSE;
                    optRotations = iterations;
                    optCutoff = cutoff;
                    initialize = 0;
                end
                
                if RMSE <= minRMSE
                    optBeta = BetaSMT;
                    optRotations = iterations;
                    optCutoff = cutoff;
                    minRMSE = RMSE;
                    lamsize = size(selLambda,1);
                    
                    
                end
            end
        end
    end
    
    X=circshift(X,[n2,0]);
    y=circshift(y,n2);
end