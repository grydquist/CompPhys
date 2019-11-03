dat = load('an0q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x,y);

figure;

dat = load('bn1q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x,y);


figure;

dat = load('an1q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x,y);

figure;
dat = load('bn2q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x,y);


figure;
dat = load('an2q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x,y);


figure;
dat = load('bn10q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x(1:end/2),y(1:end/2));


figure;
dat = load('an10q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x(1:end/2),y(1:end/2));


figure;
dat = load('an25q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x,y);


figure;
dat = load('bn25q25.txt');
x = dat(:,1);
y = dat(:,2);
plot(x,y);
