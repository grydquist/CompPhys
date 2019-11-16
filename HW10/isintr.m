function f = isintr(b)
% f = isintr(b)
% Inverse Discrete Sine Transform
% given coefficients b=[b_1,...,b_{N-1}] of sine series sum_{k=1...N-1} b_k sin(k*pi*x)
% find values f=[f_1,...,f_{N-1}] at x=1/N,...,(N-1)/N
%
% b MUST be a COLUMN VECTOR, f will also be a column vector
% If b is an array the transform acts on each column, and f is an array of the same size

N = size(b,1) + 1;
m = size(b,2);
z = zeros(1,m);
Gh = [z;b;z;-b(end:-1:1,:)]/(2i);   % use sin(t) = (exp(it)-exp(-it))/(2i)
G = ifft(Gh,2*N,1)*(2*N);           % apply ifour() to each column
f = G(2:N,:);
end