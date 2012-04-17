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

%%This is a method that uses a Minimum Description Length parameter to
%%reduce dimensionality and uses PCR to make an estimate on the smaller
%%dimension

clc;
clear;
close all;

dataset = 2;

while dataset == 1
load('Inorfull.mat');
Xraw = DATA;
Yraw = CONC;
dataset = 0;
n = size(Xraw,1); %Number of samples
X = Xraw - ones(n,1) * mean(Xraw,1);
y = Yraw - ones(n,1) * mean(Yraw,1)./(ones(n,1) * std(Yraw,0,1));
end

while dataset == 2
load('spectra.mat');
Xraw = DataSet.Spectra;
Yraw = DataSet.Conc;
dataset = 0;
n = size(Xraw,1); %Number of samples
X = Xraw - ones(n,1) * mean(Xraw,1);
y = Yraw - ones(n,1) * mean(Yraw,1);
end

while dataset == 3
load('ArtificialDataC.mat');
Xraw = Spectra;
Yraw = Conc;
dataset = 0;
n = size(Xraw,1); %Number of samples
X = Xraw - ones(n,1) * mean(Xraw,1);
y = Yraw - ones(n,1) * mean(Yraw,1);
end

while dataset == 4
load('ArtificialDataS.mat');
Xraw = Spectra;
Yraw = Conc;
dataset = 0;
n = size(Xraw,1); %Number of samples
X = Xraw - ones(n,1) * mean(Xraw,1);
y = Yraw - ones(n,1) * mean(Yraw,1);
end

while dataset == 5
load grass
%load water;
n=80;  % 50, 100, or 200
s=3;


W=pixels;

[p N]=size(W);

% compute 'true' R
avgW=mean(W,2);
W=W-avgW*ones(1,N);
Rw=1/N*W*W';
[Ew,Lambdaw] = eig(Rw);
[lambdaw,I] = sort(diag(Lambdaw),'descend');
Ew=Ew(:,I);
Lambdaw=diag(lambdaw);


%t=randn(1,p);
t=Ew(:,170)';
t=s*t/norm(t);

% Generate Gaussian samples from covariance R

W = randn([n,p])*sqrt(Lambdaw)*Ew';

%training data
sigma=1;
ytr=randn(n,1)*sigma;
Xtr=ytr*t+W;

%test data
nt=300;
yte=randn(nt,1);
Wte=randn([nt,p])*sqrt(Lambdaw)*Ew';
Xte=yte*t+Wte;
dataset = 0;

n = size(Xraw,1); %Number of samples
X = Xraw - ones(n,1) * mean(Xraw,1);
y = Yraw - ones(n,1) * mean(Yraw,1);
end

initialize = 1;

%1 - reduction every iteration
%2 - reduction at the end
reduce = 1; 

%0 No PCR
%1 Yes PCR
PCR = 1;

t = 1;
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

for iter=1:t
    iter

%      Xtr =   X(1:(t-1)*n1,:);
%      ytr =   y(1:(t-1)*n1,:);
%      Xte =   X((t-1)*n1+1:n,:);
%      yte =   y((t-1)*n1+1:n,:);

    trsize   =   size(X,1)/t;
    cal      =   1:trsize            ;
    val      =   (length(cal) + 1):size(X,1);
    var      =   1:size(X,2)       ;
    Xtr      =   X(cal,var)  ;
    ytr      =   y(cal,:)    ;
    Xte      =   X(val,var)  ;
    yte      =   y(val,:)    ;

    [ntr, mtr] = size(Xtr);
    [nte, mte] = size(Xte);
    
    MDL = 1 - exp(-(log(ntr) + 5 * log(mtr)) / ntr);
    
    NumRotations = ceil(log2(mtr)*ntr);
    
    deltaCount = 1;
    for delta = 0:0.05:1
        R=1/mtr*Xtr'*Xtr;

        T = eye(size(R));
        N = max(size(R));
        F = Criteria(R);
        [f,Irow]=max(F);
        
        iterations = 1;
        itercount = 1;  
        clear SMTArray;
        
        while iterations < NumRotations
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
                        
            if reduce == 1 || reduce == 2
                if R(i,i) < R(j,j)
                    R(i,:) = 0;
                    R(:,i) = 0;
                    choice = 1;
                else
                    R(j,:) = 0;
                    R(:,j) = 0;
                    choice = 2;
                end
            end
            
            % T=T*A
            ti = T(:,i)*cos(theta)+T(:,j)*sin(theta);   % T*A
            tj = -T(:,i)*sin(theta)+T(:,j)*cos(theta);
            T(:,i) = ti;
            T(:,j) = tj;
            
            if reduce == 1
                if choice == 1;
                   T(i,:) = 0;
                   T(:,i) = 0;
                else
                   T(j,:) = 0;
                   T(:,j) = 0;
                end
            end
            
            Fdelta = F(i,j)^2 / ((1 + delta)^2); %Correlations on the main diagonal are 1. Also, "SMT-F" has 0 for main diag
            diff(iterations,deltaCount,iter) = Fdelta - MDL; %Difference between Fdelta and MDL crit
            

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
            
            if reduce == 2
                SMTArray(iterations,1)=i;
                SMTArray(iterations,2)=j;
                SMTArray(iterations,3)=theta;
            end
            
            if Fdelta < MDL
                terminal(deltaCount) = iterations;
                iterations = NumRotations; %Breaks out of the while loop for the next rotation - this rotation is ok
            end
            
            if (Fdelta < MDL)
                [Lambda,I] = sort(lambda);
                E = T(:,I);
                Lambda = diag(Lambda);

                Lambda = diag(Lambda);                
                [Lambda,I] = sort(Lambda,'descend');
                E = E(:,I);

 
                Ered = E;
                Lambdared = Lambda;
                
                if reduce == 2
                    Ered(:,SMTArray(:,2)) = [];
                    Lambdared = Lambda(Lambda>0);

                elseif reduce == 1
                    for ecol = 1:size(E,2)
                        colsum = 0;
                        for erow = 1:size(E,1)
                            colsum = abs(E(erow,256 + 1 - ecol)) + colsum;
                        end
    
                        if colsum < 1e-10 %Account for machine (im)precision
                            Ered(:, 256 + 1 - ecol) = [];
                        end
                    end
                    Lambdared = Lambda(Lambda>0);                 
                    
                end
                
                if size(Lambdared,1) ~= size(Ered,2)
                    disp('Lambda and E are different sizes, something is wrong')
                    return
                end
                

                Lambda_inv = Lambdared.^(-1);
                
                if PCR == 0
                    
                    P = E*diag(Lambda_inv.^(1/2));
                    Ttr = Xtr*P;
                    Tte = Xte*P;
                    btilde = pinv(Ttr) * ytr;
                    
                    BetaSMT = P*btilde;
                    
                    yhat = Tte*btilde;
                    
                    RMSE = sqrt(sum(sum(((yhat - yte).^2)/(numel(yte)))));
                    
                    RMSEmat(iterations, deltaCount, iter) = RMSE;
                    
                    if initialize == 1
                        minRMSE = RMSE;
                        initialize = 0;
                    end
                    
                    if RMSE <= minRMSE
                        optBeta = BetaSMT;
                        optDelta = delta;
                        minRMSE = RMSE;
                    end
                
                elseif PCR == 1
                    Xtilde = Xtr*Ered*diag(Lambda_inv.^(1/2));
                                    
                    [U,S,P] = svd(Xtilde,'econ');
                    lambda = diag(S);
                    Tpcr = U*diag(lambda);
                    
                    pc = 20 % select number PCs (model dimension)
                    
                    for dim = 1:pc
                        dim
                        Ttr = Tpcr(:,1:dim);
                        Ptr = P(:,1:dim);
                        
                        beta = (Ttr' * Ttr) \ Ttr' * ytr;
                        betaX = Ptr * beta;
                        betaProj = Ered * diag(Lambda_inv.^(1/2)) * betaX;
                        yhat = Xte * betaProj;
                        
                        RMSE = (trace((yte - yhat)' *(yte - yhat)) / numel(yte)).^(1/2) ;
                        RMSEraw(dim, deltaCount, iter) = RMSE;
                        
                        if initialize == 1
                            minRMSE = RMSE;
                            initialize = 0;
                        end
                        
                        if RMSE <= minRMSE
                            optBeta = betaProj;
                            optDelta = delta;
                            minRMSE = RMSE;
                            optDim = dim;
                        end
                        
                    end
                end
                
            end
            
            iterations = iterations + 1;
        end
        
        deltaCount = deltaCount + 1;
    end
    

    X=circshift(X,[n2,0]);    
    y=circshift(y,n2);
end