%%This file loads data from wastewater treatment facilities EnviroX-Y.csv:
% 1-1
% 2-1
% 2-2
% 3-1
% 3-2
% 3-3
% 3-4
% 4-1
% 4-2
% 4-3
%
% It then replaces <null> entries with NaN and removes columns without
% complete input data and input rows without complete output data
%
% We then call
% Industry Standard
% PCR
% PLS
% Ridge
% OLS
%
% Cross validation takes place by choosing the most parsimonious model
% within one standard deviation of the best performing (Hastie p244)
%
% And return various quantities, most importantly optimal beta coefficients and root mean square error


clc;
clear;
close all;


%Choose the data set EnviroX-Y.csv that you want
comp1 = 3;
comp2 = 4;
filename = [ 'Enviro' num2str(comp1) '-' num2str(comp2) '.csv' ];
DATA    =   csvread(filename)  ; %Don't forget to change save at the bottom if your'e creating data!

disp('Press any key if you forgot to save variables (and you want to). Graphs are going to be generated and saved to this directory.')
pause

% Number of components

if comp1 == 1
    comp = 4;
    
elseif comp1 == 2
    if comp2 == 1
        comp = 3;
    elseif comp2 == 2
        comp = 4;
    end
    
elseif comp1 == 3
    if comp2 == 1 || comp2 == 2
        comp = 5;
    elseif comp2 == 3 || comp2 == 4
        comp = 3;
    end
    
elseif comp1 == 4
    if comp2 == 1 || comp2 == 2
        comp = 7;
    elseif comp2 == 3
        comp = 4;
    end
end

Spectra     =  DATA(2:end, 2:(end - comp))   ;
Conc        =  DATA(2:end,(end - comp + 1):end) ;

%Replace null observations with NaN
Conc(find(Conc==0)) =   NaN;

Xraw = Spectra;
Yraw = Conc;

%Number of cross validation blocks
t = 3;

%Find where we don't have input data
xmatindex = ~isnan(Xraw);
for col = 1:size(Xraw,2)
    if sum(xmatindex(:,col)) == size(Xraw,1)
        xindex(1,col) = 1;
    else
        xindex(1,col) = 0;
    end
end

%Remove columns without input data
Xnan = Xraw(:, xindex == 1);
outputs = size(Yraw,2);

%One component at a time
for i = 1:outputs 
    i
    ycol = Yraw(:,i);
    
    %Remove input/output data where we don't have output data
    yindex = ~isnan(ycol);
    X = Xnan(yindex == 1, :);
    nx = size(X,1);
    X = (X - ones(nx,1) * mean(X)) / std2(X);
    ynan = ycol(yindex == 1,:);
    
    y = (ynan - mean(ynan)) / std(ynan);
    
    [y, yindex] = sort(y,'descend');
    X = X(yindex,:);
    
    
    sizex = size(X);
    sizey = size(y);
    
    
    %Industry standard
    [optBeta, RMSEmat, stddev, optdim, rsqall] = industry(X, y, t);
    
    RMSEmat = (RMSEmat) * std(ynan);
    
    %Size of the coefficients depends on the size of the outputs, so you
    %have to insert NaNs to make the matrix the right size
    if i ~= 1 && size(optBeta,1) < size(betaIndu,1)
        optBeta((size(optBeta,1) + 1):size(betaIndu,1), :) = NaN;
    end
    if i ~= 1 && size(optBeta,1) > size(betaIndu,1)
        betaIndu((size(betaIndu,1) + 1):size(optBeta,1), :) = NaN;
    end
    
    betaIndu(:, i) = optBeta;
    
    if i ~= 1 && size(RMSEmat,1) < size(rmseIndu,1)
        RMSEmat((size(RMSEmat,1) + 1):size(rmseIndu,1), :) = NaN;
    end
    if i ~= 1 && size(RMSEmat,1) > size(rmseIndu,1)
        rmseIndu((size(rmseIndu,1) + 1):size(RMSEmat,1), :) = NaN;
    end
    
    rmseIndu(:, (3*i - 2):(3*i - 3 + size(RMSEmat,2))) = RMSEmat;
    dimindu(i) = optdim;
    
    if i ~= 1 && size(stdIndu,1) < size(stddev,1)
        stdIndu((size(stdIndu,1) + 1):size(stddev,1), :) = NaN;
    end
    if i ~= 1 && size(stdIndu,1) > size(stddev,1)
        stddev((size(stddev,1) + 1):size(stdIndu,1), :) = NaN;
    end
    
    stdIndu(:,i) = stddev * abs(std(ynan));
    
    
    %Size of the coefficients depends on the size of the outputs, so you
    %have to insert NaNs to make the matrix the right size
    %PCR
    [optBeta, RMSEmat, optdim, R2, stddev] = PCR(X, y, t);
    betaPCR(:, i) = optBeta;
    yPCR(1:size(X,1), i) = X * optBeta;
    RMSEmat = (RMSEmat) * std(ynan);
    
    if i ~= 1 && size(RMSEmat,1) < size(rmsePCR,1)
        RMSEmat((size(RMSEmat,1) + 1):size(rmsePCR,1), :) = NaN;
    end
    if i ~= 1 && size(RMSEmat,1) > size(rmsePCR,1)
        rmsePCR((size(rmsePCR,1) + 1):size(RMSEmat,1), :) = NaN;
    end
    rmsePCR(:, (3*i - 2):(3*i - 3 + size(RMSEmat,2))) = RMSEmat;
    dimPCR(i) = optdim;
    if i ~= 1 && size(R2,1) < size(pcrr2,1)
        R2((size(R2,1) + 1):size(pcrr2,1), :) = NaN;
    end
    if i ~= 1 && size(R2,1) > size(pcrr2,1)
        pcrr2((size(pcrr2,1) + 1):size(R2,1), :) = NaN;
    end
    
    if i ~= 1 && size(pcrstd,1) < size(stddev,1)
        pcrstd((size(pcrstd,1) + 1):size(stddev,1), :) = NaN;
    end
    if i ~= 1 && size(pcrstd,1) > size(stddev,1)
        stddev((size(stddev,1) + 1):size(pcrstd,1), :) = NaN;
    end
    pcrstd(:,i) = stddev * abs(std(ynan));
    pcrr2(:, i) = R2;
    
    
    %PLS
    [optBeta, RMSEmat, optdim, R2, stddev] = PLS(X, y, t);
    betaPLS(:, i) = optBeta;
    yPLS(1:size(X,1), i) = X * optBeta;
    RMSEmat = (RMSEmat) * std(ynan);
    rmsePLS(:, (3*i - 2):(3*i - 3 + size(RMSEmat,2))) = RMSEmat;
    dimPLS(i) = optdim;
    plsr2(:, i) = R2;
    plsstd(:,i) = stddev * abs(std(ynan));
    
    
    %Ridge
    [optBeta, RMSEmat, optlambda, edof, minedof, R2, stddev] = Ridge(X, y, t);
    betaRidge(:, i) = optBeta;
    yRidge(1:size(X,1), i) = X * optBeta;
    edofRidge(:, (3*i - 2):(3*i - 3 + size(edof,2))) = edof;
    minEdof(:, i) = minedof;
    RMSEmat = (RMSEmat) * std(ynan);
    rmseRidge(:, (3*i - 2):(3*i - 3 + size(RMSEmat,2))) = RMSEmat;
    lambdaRidge(i) = optlambda;
    ridgestd(:,i) = stddev * abs(std(ynan));
    ridger2(:,i) = R2;
end

%Make graphs!
rmsrerror(rmsePCR, pcrstd, rmsePLS, plsstd, rmseRidge, ridgestd, rmseIndu, stdIndu, comp1, comp2, outputs)

% Uncomment if you want to save data sets
% 
% savefile = ['PCR' num2str(comp1) '-' num2str(comp2) '.mat' ];
% save(savefile, 'betaPCR' , 'yPCR' , 'rmsePCR' , 'dimPCR' , 'pcrr2')
% 
% savefile = ['PLS' num2str(comp1) '-' num2str(comp2) '.mat' ];
% save(savefile, 'betaPLS' , 'yPLS' , 'rmsePLS' , 'dimPLS' , 'plsr2')
% 
% 
% savefile = ['Ridge' num2str(comp1) '-' num2str(comp2) '.mat' ];
% save(savefile, 'betaRidge' , 'yRidge' , 'edofRidge' , 'minEdof' , 'rmseRidge' , 'lambdaRidge' , 'ridger2')
% 
% 
% savefile = ['industry' num2str(comp1) '-' num2str(comp2) '.mat' ];
% save(savefile, 'betaIndu' , 'rmseIndu' , 'dimindu')