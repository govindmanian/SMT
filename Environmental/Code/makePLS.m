%%PLS1 algorithm for PLS

function model = makePLS(X,Y,kmax)  

nx  =   size(X,2)   ;
A   =   X'*Y        ;
M   =   X'*X        ;
C   =   eye(nx)     ;

for k = 1:kmax
    [u,s,v] =   svds(A'*A,1)    ;
    q       =   u(:,1)          ;
    w       =   C*A*q           ; 
    w       =   w/norm(w)       ;
    W(:,k)  =   w               ;
    p       =   M*w             ;
    c       =   w'*M*w          ;
    p       =   p/c             ;
    P(:,k)  =   p               ;
    q       =   A'*w/c          ;
    Q(:,k)  =   q               ;
    A       =   A - c*p*q'      ;
    M       =   M - c*p*p'      ;
    C       =   C - w*p'        ;
end

model.type  =   'PLS'   ;
model.P     =   P       ;
model.W     =   W       ;
model.Q     =   Q       ;   