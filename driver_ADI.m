% Deriver code for 
% ADI for solving matrix equation AX - XB = C for small matrix
% and compare its performance with fADI

%% 
% initialize matrix A, B as normal matrices
n = 100;
a = 1; b = 4;  % interval for eigs of A is [a, b]
[Q, ~] = qr(rand(n,n));
Aeigs = (b-a)*rand(n,1) + a; 
A = Q * diag(Aeigs) * Q';
B = -A; 
Beigs = diag(B);% B is reflection of A across imaginary axis

% initialize right-hand-side matrix C = GF' of rank rnk_C
rnk_C = 3;
G = rand(n,rnk_C); F = rand(n,rnk_C); 
C = G*F';

%%
% get shift paramters
% for k steps of ADI or fADI
k = 12; 
I = [[a,b], [-b, -a]];
[p, q] = getshifts_adi(I,k); 

% use Bartels-Stewart to solve AX - XB = C
% seen as 'precise' solution
X_ly = lyap(A, -B, -C);

%% call ADI and check relative 2 norm error and norm of residual
X_ADI = adi(A, B, C, p, q);
rel2err_ADI = norm(X_ly - X_ADI)/norm(X_ly)
residual_ADI = norm(A*X_ADI - X_ADI*B - C)

%% call fADI and check relative 2 norm error and norm of residual
[ZZ, DD, YY] = fadi(A,B,G,F, p,q);
X_fADI = ZZ*DD*YY';
rel2err_fADI = norm(X_ly - X_fADI)/norm(X_ly)
residual_fADI = norm(A*X_fADI - X_fADI*B - C)