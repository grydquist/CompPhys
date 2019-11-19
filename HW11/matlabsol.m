fi = @(x) exp(-x.^2);
fa = @(x,t) 1/sqrt(1+4*t).*exp(-x.^2./(1+4*t));
L = 10;
nx = 100;
dt = 0.001;
Nt = 1000;
x = linspace(0,L,nx);
dx = x(2)-x(1);
t=0;
u = fb(x);
ut = u;
ua = fa(x,t);
plot(x,u,'o',x,ua);

%% FTCS

% for i = 1:Nt
%     t = t+dt;
%     for j=1:nx
%         if j==1
%             ut(j) = u(j) +dt*(2*u(j+1)-2*u(j))/dx^2;
%         elseif j==nx
%             ut(j) = u(j) +dt*(u(j-1)-2*u(j))/dx^2;
%         else
%             ut(j) = u(j) +dt*(u(j+1)-2*u(j)+u(j-1))/dx^2;
%         end
%     end
%     u = ut;
%     ua = fa(x,t);
%     clf;
%     plot(x,u,'o',x,ua);
%     axis([0,L,0,1]);
%     pause(0.01);
% end

%% CN
A = zeros(nx,nx);
dt = 0.1;
A(1:nx+1:end) = (1/dt+1/dx^2);
A(2:nx+1:end) = -1/2/dx^2;
A = A';
A(2:nx+1:end) = -1/2/dx^2;
A(1,2) = 2*A(1,2);

b = zeros(nx,1);

for i = 1:Nt
    t = t+dt;
    for j=1:nx
        if j==1
            b(j) = u(j)/dt+(2*u(j+1)-2*u(j))/2/dx^2;
        elseif j==nx
            b(j) = u(j)/dt+(-2*u(j)+u(j-1))/2/dx^2;
        else
            b(j) = u(j)/dt+(u(j+1)-2*u(j)+u(j-1))/2/dx^2;
        end
    end
    u = A\b;
    ua = fa(x,t);
    clf;
    plot(x,u,'o',x,ua);
    axis([0,L,0,1]);
    pause(0.01);
end

%%
d = load('CN.txt');
u = reshape(d(:,1),100,1000);
ua = reshape(d(:,2),100,1000);
x  = reshape(d(:,3),100,1000);
x = x(:,1);
for i = 1:Nt
    clf;
    plot(x,ua(:,i),'o',x,ua(:,i));
    axis([0,L,0,1]);
    pause(0.01);
    
end
