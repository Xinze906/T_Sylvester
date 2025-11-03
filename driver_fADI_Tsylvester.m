% Driver code for solving AX+X'B = C by fADI 
% via (B^{-T}A) X - X(A^{-T} B) = B^{-T} C - B^{-T} C^{T} A{-T} B

%%
% initialize A and B as normal matrices
n = 100;
Aeigs = linspace(10, 12, n);
Beigs = 50 + rand(n, 1)*10;
[Q, ~] = qr(rand(n,n));

A = Q * diag(Aeigs) *Q';
B = Q * diag(Beigs) *Q';

% initialize C as a random matrix of rank rnk_C
rnk_C = 3;
C1 = rand(n, rnk_C);
C2 = rand(n, rnk_C);
C = C1 * C2';

%%
% seperate spec(B^{-T}A) and spec(A^{-T} B) using shift parameters on disks
% must be disjoint
% specBA = eig(B' \ A);
% specAB = eig(A' \ B);
% 
% oBA = (specBA(1, 1) + specBA(end, 1))/2;
% rBA = max(abs(specBA(end, 1) - oBA), abs(specBA(1, 1) - oBA));
% 
% oAB = (specAB(1, 1) + specAB(end, 1))/2;
% rAB = max(abs(specAB(end, 1) - oAB), abs(specAB(1, 1) - oAB));
% 
% [p, q] = getshifts_smith([oBA, rBA, oAB, rAB]);
%%
% get shift parameters
% seperate spec(B^{-T}A) and spec(A^{-T} B) using shift parameters on intervals
% must be disjoint

specBA = eig(B' \ A);
specAB = eig(A' \ B);

a = min(specBA); b = max(specBA);
c = min(specAB); d = max(specAB);
I = [a b c d];
num_shift = 12;
[p, q] = getshifts_adi(I, num_shift);
% prepare T-Sylvester in the form of Sylvester
U = [B'\C1, B'\C2];
V = [C2'; -C1'*(A'\B)];

%%
% solve via fADI, check residual and relative 2 norm error
[ZZ, DD, YY] = fadi(B' \ A, A' \ B,  U, V', p, q);
X_approx = ZZ * DD * YY';

% res for T-Sylvester
res_2_2 = norm(A * X_approx + X_approx' * B - C); 
% res for the form we solved
res_2_3 = norm((B' \ A) * X_approx - X_approx * (A' \ B) - U*V); 
% res for another form
res_2_4 = norm(A * X_approx * A' - B' * X_approx * B - B' * U*V);

% relative 2 norm error
X_real = lyap( B'\A, -A'\B, -U*V);
rel2err = norm(X_real-X_approx) / norm(X_real);