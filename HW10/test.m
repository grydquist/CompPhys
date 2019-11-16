A=1; B=1;                      % lengths of rectangle sides
fct = @(x) x*0+1;     % given function for top side of rectangle y=B
N = 16;                        % use spacing A/N for x
x = A*(1:N-1)'/N;              % only interior nodes            (column vector)
g = fct(x);                    % evaluate g for 1/N,...,(N-1)/N (column vector)


b = sintr(g);                  % take Discrete Sine Transform
kv = (1:N-1)';                 % column vector of k-values
muv = kv*pi/A;                 % column vector of mu-values
C = b./sinh(muv*B);            % coefficients C_k from solution formula

M = 30;
y = B*(0:M)/M;                 % find solution for these y-values (row vector)

W = sinh(muv*y);               % array of sinh(mu_k*y) for k=1...N-1, for all y-values
                               % (column vector)*(row vector) gives array
U = isintr(diag(C)*W);         % multiply 1st row by C(1), 2nd row by C(2), ...
                               % then take Inverse Sine Transform of each column

surf(x,y,U');                  % make surface plot of solution
xlabel('x'); ylabel('y'); title('u(x,y)')
view(-15,30);

function b = sintr(f)
% b = sintr(f)
% Discrete Sine Transform
% given values f=[f_1,...,f_{N-1}] at x=1/N,...,(N-1)/N (NOTE: endpoints are NOT included)
% find coefficients b=[b_1,...,b_{N-1}] of sine series sum_{k=1...N-1} b_k sin(k*pi*x)
%
% f MUST be a COLUMN VECTOR, b will also be a column vector
% If f is an array the transform acts on each column, and b is an array of the same size

N = size(f,1) + 1;
m = size(f,2);
z = zeros(1,m);
G = [z;f;z;-f(end:-1:1,:)];        % odd extension of each column of length n=N-1 to length 2*N
Gh = fft(G,2*N,1)/(2*N);           % apply four() to each column
b = 1i*(Gh(2:N,:)-Gh(end:-1:N+2,:));
end

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