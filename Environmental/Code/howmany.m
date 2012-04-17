DATA    =   csvread('Enviro4-3.csv')  ;

comp = 4;

Spectra     =   DATA(2:end, 2:(end - comp))   ; %1: 4; 4-2: 7;
Conc        =   DATA(2:end,(end - comp + 1):end) ; %1 - 2:end,219:222; 2 - 2:end,(end-3):end

Conc(find(Conc==0)) =   NaN         ;

Yraw = Conc;  
yindex = ~isnan(Yraw);

sum(yindex)