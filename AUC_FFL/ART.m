function [ X, k ] = ART(A0, b0, X0, e0)
% 重建方法
%ART:to solve linear equation Ax = b
%Input  -   A: the coefficient matrix (matrix's size is nxn)
%       -   b: the constant term(n-dimensions vector)
%       -  X0: the inital point of iteration
%       -  e0: the termination condition of iteration
%Output -   X: the solution of the linear equation(n-dimensions vector)
%       -   k: the times of itation

A  = A0;
b = b0;
X = X0;
k = 0;
e = 2; % 迭代误差
[m,n] = size(A);

while (norm(e)>e0 && k <= 100)
    for i = 1:n
        unitLen = norm(A(i,:));
        d = (b(i)-A(i,:)*X)/unitLen;
        Xf = X;
        X = X+(d.*(A(i,:)./unitLen))';
        e = norm(Xf-X);
    end
    k = k+1;
end
end