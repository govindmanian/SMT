clc
clear all
close all

DATA    =   csvread('Enviro2.csv')  ;

Sample      =   DATA(2:end,1)       ;
Wavelength  =   DATA(1,3:218)       ;
Spectra     =   DATA(2:end,3:218)   ;
Conc        =   DATA(2:end,219:222) ;

Conc(find(Conc==0)) =   NaN         ;

figure
plot(Sample,Conc)

figure
plot(Wavelength,Spectra')
