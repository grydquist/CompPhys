fileid = fopen('coeffs.txt');
c = fscanf(fileid, '%f');
fclose(fileid);

x = 0:(length(c)-1);
semilogy(x,abs(c),'o');
xlabel('i');
ylabel('|c_i|');

figure;
plot(x,c,'o')
xlabel('i');
ylabel('c_i');

fileid = fopen('derivs.txt');
a = fscanf(fileid, '%f');
fclose(fileid);

dera = a(2:2:end);
derc = a(1:2:end);
xx = linspace(-0.5,2,length(dera));

figure;
plot(xx,derc,'o');
hold on
plot(xx,dera);
xlabel('x');
ylabel('df(x)/dx');
legend('Chebyshev Approximation','Analytical','Location','southeast')