% driver file for solving AX+X'B = C by fadi method
% using (B^{-T}A) X - X(A^{-T} B) = B^{-T} C - B^{-T} C^{T} A{-T} B
clear all 
close all

% initialize A and B as normal matrices
n = 100;
Aeigs = linspace(10, 12, n);
Beigs = 50 + rand(n, 1)*10;
[Q, ~] = qr(rand(n,n));

A = Q * diag(Aeigs) *Q';
B = Q * diag(Beigs) *Q';

% initialize C as a random matrix
rnk = 1;
C1 = rand(n, rnk);
C2 = rand(n, rnk);
C = C1 * C2';

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

% seperate spec(B^{-T}A) and spec(A^{-T} B) using shift parameters on intervals
% must be disjoint


specBA = eig(B' \ A);
specAB = eig(A' \ B);


a = min(specBA); b = max(specBA);
c = min(specAB); d = max(specAB);
I = [a b c d];
num_shift = 12;
[p, q] = getshifts_adi(I, num_shift);
U = [B'\C1, B'\C2];
V = [C2'; -C1'*(A'\B)];
[ZZ, DD, YY] = fadi(B' \ A, A' \ B,  U, V', p, q);


% norm(X-ZZ*DD*YY') B-S

X_approx = ZZ * DD * YY';
residual_2_3 = norm((B' \ A) * X_approx - X_approx' * (A' \ B) - U*V);

residual_2_2 = norm(A* X_approx + X_approx' *B - C);
X_real = lyap( B'\A, -A'\B, -U*V);
error = norm(X_real-X_approx) / norm(X_real);