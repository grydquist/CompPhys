fileid = fopen('coeffs.txt');
c = fscanf(fileid, '%f');
fclose(fileid);

x = 0:(length(c)-1);
semilogy(x,abs(c),'o');

fileid = fopen('derivs.txt');
a = fscanf(fileid, '%f');
fclose(fileid);

dera = a(2:2:end);
derc = a(1:2:end);
xx = linspace(-0.5,2,length(dera));


plot(xx,derc,'o');
hold on
plot(xx,dera);
