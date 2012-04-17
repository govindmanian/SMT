clc;
clear;
close all;

DATA    =   csvread('Enviro1-1.csv')  ; %Don't forget to change save at the bottom!
load EnviroRidge1-1
load EnviroPCR1-1
load EnviroPLS1-1
load EnviroOLS1-1

% Number of components
% 1.1 4
% 2.1 3
% 2.2 4
% 3.1 5
% 3.2 5
% 3.3 3
% 3.4 3
% 4.1 7
% 4.2 7
% 4.3 4

comp = 4;
Spectra     =   DATA(2:end, 2:(end - comp))   ;
Conc        =   DATA(2:end,(end - comp + 1):end) ;

Conc(find(Conc==0)) =   NaN         ;

Xraw = Spectra;
Yraw = Conc;  

xmatindex = ~isnan(Xraw);

for col = 1:size(Xraw,2)
    if sum(xmatindex(:,col)) == size(Xraw,1)
        xindex(1,col) = 1;
    else
        xindex(1,col) = 0;
    end
end

Xx = Xraw(:, xindex == 1);

n = size(Xraw,1); %Number of samples
Xraw = (Xraw - ones(n,1) * mean(Xraw)) / std2(Xraw); %Group scaling

for i = 1:size(Yraw,2) %One component at a time
    ycol = Yraw(:,i);
    yPCRi = yPCR(:,i);
    yPLSi = yPLS(:,i);
    yRidgei = yRidge(:,i);
    yOLSi = yOLS(:,i);
    
    yPCRi(find(yPCRi==0)) = [];
    yPLSi(find(yPLSi==0)) = [];
    yRidgei(find(yRidgei==0)) = [];
    yOLSi(find(yOLSi==0)) = [];

    
    yindex = ~isnan(ycol);
    X = Xx(yindex == 1, :);
    ynan = ycol(yindex == 1,:);
        
    
    y = (ynan - mean(ynan)) / std(ynan);    
    [y, reindex] = sort(y,'descend');
    X = X(reindex,:);
    
    pPCR(i) =  corr(y, yPCRi)    ;    
    pPLS(i) = corr(y, yPLSi)    ;
    pRidge(i) = corr(y, yRidgei) ;
    pOLS(i) = corr(y, yOLSi)    ;
    

end

save Pearson1-1 pPCR pPLS pRidge pOLS
