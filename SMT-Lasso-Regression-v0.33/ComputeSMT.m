function   [T,lambda,SMTArray,likelihood] = ComputeSMT(R,NumRotations,Y)
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
% likelihood: likelihood of the test data

if nargin == 3
    flag_design = 1; % design process
    M=size(Y,2);
    S=1/M*Y*Y';
else
    flag_design =0;  % generate the desired SMT covariance estimate
end
likelihood =zeros(NumRotations,1);

T = eye(size(R));
N = max(size(R));


SMTArray = zeros(NumRotations,3);

F = Criteria(R);
[f,Irow]=max(F);
for iterations=1:NumRotations

%%%% slow search %%%%%%%
%  [f,Irow]=max(F);
%  [val,Icol]=max(f);
%  i=Irow(Icol);
%  j=Icol;

%%%% fast search %%%%%%%
  [val,Icol]=max(f);
  i=Irow(Icol);
  j=Icol;

  
  if(i>j) 
    ii = i; i=j; j=ii;
  end

  %disp(sprintf('i=%d,j=%d \n',i,j));

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
  for n=1:N
      if Irow(n)==i || Irow(n)==j || n==i || n==j
          [f(n),Irow(n)]=max(F(:,n));
      else
          if f(n)<F(i,n)
              f(n)=F(i,n);
              Irow(n)=i;
          end
          if f(n)<F(j,n)
              f(n)=F(j,n);
              Irow(n)=j;
          end
      end
  end
 

%  MDL(iterations)=M*N/2+M/2*sum(log(diag(R)))+N*M/2*log(2*pi)+iterations/2*log(N*M);

  if flag_design==1 % && mod(iterations,ceil(NumRotations/4))==0
        % S=A'*S*A 
        ci = S(:,i)*cos(theta)+S(:,j)*sin(theta);   % R*A
        cj = -S(:,i)*sin(theta)+S(:,j)*cos(theta);
        S(:,i) = ci;
        S(:,j) = cj;

        ri = S(i,:)*cos(theta) + S(j,:)*sin(theta); % A'*R*A
        rj = -S(i,:)*sin(theta) + S(j,:)*cos(theta); 
        S(i,:) = ri;  
        S(j,:) = rj;
        
        likelihood(iterations) = -1/2*sum(diag(S)./diag(R))-1/2*sum(log(diag(R)))-N*1/2*log(2*pi);
  end
          
  
  %%%%% sparse representation of Givens rotations %%%%%%
  SMTArray(iterations,1)=i;
  SMTArray(iterations,2)=j;
  SMTArray(iterations,3)=theta;
  lambda=diag(R);
end


