function blasso=LassoRegression(y,X,method,option);

if strcmp(method,'CV')
    [n,m]=size(y);
    %n is the number of samples;
    %m is the number of "outputs"
    [n,p]=size(X);

    blasso=zeros(p,m);

    for i=1:m;
        blasso(:,i)=LassoCV(y(:,i),X,method,option);
    end
else
   fprintf('Wrong method option\n');
   bshrk=0;
   value=0;
end