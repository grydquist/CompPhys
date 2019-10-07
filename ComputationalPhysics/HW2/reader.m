x5=zeros(1,300);
y5=x5;
x10 = x5;
y10 = x5;
x20 = x5;
y20 = x5;
full5 = load('p5.txt');
full10 = load('p10.txt');
full20 = load('p20.txt');
for i =1:300
    x5(i) = full5(2*i-1);
    y5(i) = full5(2*i);
    x10(i) = full10(2*i-1);
    y10(i) = full10(2*i);
    x20(i) = full20(2*i-1);
    y20(i) = full20(2*i);
end
ya = 1./(1+25*x5.^2);
plot(x5,ya,x5,y5,'--r')
title('5 Points (polynomial)')
legend('Analytical','Polynomial')
figure
plot(x5,ya,x10,y10,'--r')
title('10 Points (polynomial)')
legend('Analytical','Polynomial')
figure
plot(x5,ya,x20,y20,'--r')
title('20 Points (polynomial)')
legend('Analytical','Polynomial')
figure

EN=zeros(1,3);
EN(1) = max(abs(ya-y5));
EN(2) = max(abs(ya-y10));
EN(3) = max(abs(ya-y20));
plot([5,10,20],EN)
xlabel('N')
ylabel('E_N')
title('Error vs. Polynomial Degree')
figure


full5 = load('c5.txt');
full10 = load('c10.txt');
full20 = load('c20.txt');
for i =1:300
    x5(i) = full5(2*i-1);
    y5(i) = full5(2*i);
    x10(i) = full10(2*i-1);
    y10(i) = full10(2*i);
    x20(i) = full20(2*i-1);
    y20(i) = full20(2*i);
end

plot(x5,ya,x5,y5,'--r')
title('5 Points (cubic spline)')
legend('Analytical','Cubic spline')
figure
plot(x5,ya,x10,y10,'--r')
title('10 Points (cubic spline)')
legend('Analytical','Cubic spline')
figure
plot(x5,ya,x20,y20,'--r')
title('20 Points (cubic spline)')
legend('Analytical','Cubic spline')
figure

EN=zeros(1,3);
EN(1) = max(abs(ya-y5));
EN(2) = max(abs(ya-y10));
EN(3) = max(abs(ya-y20));
plot([5,10,20],EN)
xlabel('N')
ylabel('E_N')
title('Error vs. Number of points')