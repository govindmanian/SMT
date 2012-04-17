function F = criteria(R)

lambda = diag(R);

F = (R-diag(lambda)).^2 ./ (lambda*lambda');

