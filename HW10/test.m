A=1; B=1;                      % lengths of rectangle sides
fct = @(x) x*0+1;     % given function for top side of rectangle y=B
N = 8;                        % use spacing A/N for x
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

%surf(x,y,U');                  % make surface plot of solution
%xlabel('x'); ylabel('y'); title('u(x,y)')
%view(-15,30);

rho = zeros(N,N);
rho(1,:) = -1;
aa = sintr(rho');
rhoh = sintr(aa');


%rhoh = dstn(rho);

m= 1:N;
n = 1:N;
[NN,MM] = meshgrid(n,m);

uh = rhoh./2./(cos(pi*MM/N)+cos(pi*NN/N)-2);
aa = sintr(uh')*N/2;
u  = sintr(aa')*N/2;

surf(U)

