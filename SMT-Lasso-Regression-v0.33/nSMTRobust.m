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

%%This is a version of SMT that provides robust covariance estimates by
%%using an adaptation of Ridge Regularization

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


t = 4;
n1 = ceil(n*1/t);
n2 = n - (t-1)*n1;

totRMSE = 0; %Numerator for average RMSE
totM = 0;
avgRMSE = 0;

for iter=1:t

    trsize   =   100;
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
    RMSEcount = 1;

   
    R=1/mtr*(Xtr')*Xtr;

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
        
        if mod(iterations, NumRotations) == 0
            [Lambda,I] = sort(lambda);
            E = T(:,I);
            Lambda = diag(Lambda);

            Lambda = diag(Lambda);                
            [Lambda,I] = sort(Lambda,'descend');
            E = E(:,I);

            %Selection
            sel = 256;
            E = E(:, 1:sel);
            Lambda = Lambda(1:sel);

            Lambda_inv = Lambda.^(-1);
            P = E*diag(Lambda_inv.^(1/2));

            Rsmt    = pinv(P)'* diag(ones(sel,1)*mtr) *pinv(P);
            BetaSMTridge = (Rsmt) \ Xtr' * ytr;

            yhat = Xte * BetaSMTridge;
            RMSE = sqrt(sum(sum(((yhat - yte).^2)/(numel(yte)))));
            
            RMSEmat(iter) = RMSE;
            
            if (iter == 1) && (RMSEcount == 1)
                minRMSE = RMSE;
                optRotations = iterations;
            end

            if RMSE <= minRMSE
                optBeta = BetaSMTridge;
                optRotations = iterations;            
                minRMSE = RMSE;
            end
            
        end

    end

        X=circshift(X,[n2,0]);
        y=circshift(y,n2);
end


%%%Now take the best beta coefficients and see how you did!

% finyhat = finXte * optBeta
%finish coding

makefig = 0;
while makefig == 1
    figure
    for iter=1:t
        randomColor = rand(1,3);
        plot(1:NumRotations, RMSEmat(:,iter+1), '.', 'Color', randomColor);
        title('Number of Rotations vs RMSE');
        xlabel('Number of Rotations');
        ylabel('RMSE');
        if iter == 1
            hold on;
        end
    end
    hold on
    avgmat(1:NumRotations,1) = avgRMSE;
    plot(1:NumRotations, avgmat, '--r')
    title('RMSE for each dimension');
    ylabel('RMSE');
    xlabel('Dimension');
    hold off

    figure
    index = [1:1:size(optBeta,1)]';
    plot(index, optBeta(:,1), 'o', index, optBeta(:,2), 'o');
    title('Optimal beta coefficients');
    ylabel('Magnitude');
    xlabel('Index');
    makefig = 0;
end

