function [a,b,theta] = rotation(rii,rjj,rij)

tmp1 = 0.5*(rii + rjj);
tmp2 = 0.5*sqrt( (rii-rjj)^2 + 4*rij^2 );

a = tmp1 + tmp2;
b = tmp1 - tmp2;

theta = (atan2(2*rij,rii-rjj))/2;

