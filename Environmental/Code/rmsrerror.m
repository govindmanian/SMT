function [] = rmsrerror(rmsePCR, pcrstd, rmsePLS, plsstd, rmseRidge, ridgestd, rmseIndu, stdIndu, comp1, comp2, outputs)



for i = 1:outputs
    i
    l = 3*i - 2;
    
    %%% PCR
    
    filename = [ 'rmse-pcr-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.fig' ];
    
    [pcrmin, index] = min(rmsePCR(:,l));
    pcrmax = pcrmin + pcrstd(index,i);
    
    donezo = 0;
    dimpcr = 1;
    for row = 1:index
        if rmsePCR(row,l) < pcrmax && donezo == 0
            dimpcr = row;
            donezo = 1;
        end
    end
    
    if sum(isnan(rmsePCR(:,l))) >= 1
        notnan = ~isnan(rmsePCR(:,l));
        pcrtemp = rmsePCR(notnan,l);
        notnan = ~isnan(pcrstd(:,i));
        E = pcrstd(notnan,i);
        h = errorbar(pcrtemp, E);
    else
        E = pcrstd;
        h = errorbar(rmsePCR(:,l), pcrstd(:,i));
    end
    hold on
    plot(dimpcr, rmsePCR(dimpcr,l), 'r*');
    title(filename)
    xlabel('Number of dimensions')
    ylabel('Magnitude of Root Mean Square Error')
    saveas(h,filename)
    hold off
    
    
    filename = [ 'rmse-pcr-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.mat' ];
    gpcr = [dimpcr, rmsePCR(dimpcr,l)];
    save(filename, 'gpcr')
    
    %%% PLS
    [plsmin, index] = min(rmsePLS(:,l));
    plsmax = plsmin + plsstd(index,i);
    
    donezo = 0;
    dimpls = 1;
    for row = 1:index
        if rmsePLS(row,l) < plsmax && donezo == 0
            dimpls = row;
            donezo = 1;
        end
    end
    
    filename = [ 'rmse-pls-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.fig' ];
    plot(dimpls, rmsePLS(dimpls,l), 'r*');
    hold on
    h = errorbar(rmsePLS(:,l), plsstd(:,i));
    title(filename)    
    xlabel('Number of dimensions')
    ylabel('Magnitude of Root Mean Square Error')
    saveas(h,filename)
    hold off
    
    filename = [ 'rmse-pls-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.mat' ];
    gpls = [dimpls, rmsePLS(dimpls,l)];
    save(filename, 'gpls')
    
    
    %%% INDUSTRY
    
    filename = [ 'rmse-industry-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.fig' ];
    
    [indumin, index] = min(rmseIndu(:,l));
    indumax = indumin + stdIndu(index,i);
    
    
    donezo = 0;
    dimindu = 1;
    for row = 1:index
        if rmseIndu(row,l) < indumax && donezo == 0
            dimindu = row;
            donezo = 1;
        end
    end
    
    if sum(isnan(rmseIndu(:,l))) >= 1
        notnan = ~isnan(rmseIndu(:,l));
        indutemp = rmseIndu(notnan,l);
        notnan = ~isnan(stdIndu(:,i));
        E = stdIndu(notnan,i);
        h = errorbar(indutemp, E);
    else
        E = stdIndu;
        h = errorbar(rmseIndu(:,l), stdIndu(:,i));
    end
    hold on    
    plot(dimindu, rmseIndu(dimindu,l), 'r*');
    title(filename)    
    xlabel('Number of dimensions')
    ylabel('Magnitude of Root Mean Square Error')
    saveas(h,filename)
    hold off
    
    filename = [ 'rmse-industry-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.mat' ];
    gindu = [dimindu, rmseIndu(dimindu,l)];
    save(filename, 'gindu')
    
    %%% RIDGE
    
    [ridgemin, index] = min(rmseRidge(:,l));
    ridgemax = ridgemin + ridgestd(index,i);
    
    donezo = 0;
    dimridge = 1;
    for row = 1:index
        if rmseRidge(row,l) < ridgemax && donezo == 0
            dimridge = row;
            donezo = 1;
        end
    end
    
    filename = [ 'rmse-ridge-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.fig' ];
    plot(dimridge, rmseRidge(dimridge,l), 'r*');
    hold on
    h = errorbar(rmseRidge(:,l), ridgestd(:,i));
    title(filename)    
    xlabel('Magnitude of Penalty')
    ylabel('Magnitude of Root Mean Square Error')
    saveas(h,filename)
    hold off
    
    filename = [ 'rmse-ridge-' num2str(comp1) '-' num2str(comp2) '-' num2str(i) '.mat' ];
    gridge = [dimridge, rmseRidge(dimridge,l)];
    save(filename, 'gridge')
end


end