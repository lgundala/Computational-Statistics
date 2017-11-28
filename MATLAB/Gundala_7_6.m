clc;
clear;
%%
%gibbs sampling
A = importdata('coal.dat');
y = A.data(:,2);
n = length(y);
m = 1000;
a1 = 3;
a2 = 3;
lambda1 = zeros(1,m);
lambda2 = zeros(1,m);
theta = zeros(1,n);
like = zeros(1,n);
theta(1) = unifrnd(0,100);
lambda1(1) = 1;
lambda2(1) = 1;
b1 = 1;
b2 = 1;
c1 = 0;
c2 =0;
d1 = 1;
d2 = 1;
for i = 2:m
    kk = theta(i-1);
    t = a1+sum(y(1:kk));
    lam = theta(i-1);
    lambda1(i) = gamrnd(t,1/lam,1,1);
    t = a2+sum(y) - sum(y(1:kk));
    lam = n-kk+b2;
    lambda2(i) = gamrnd(t,1/lam,1,1);
    b1 = gamrnd(a1+c1,1/(lambda1(i)+d1),1,1);
    b2 = gamrnd(a2+c2,1/(lambda1(i)+d2),1,1);
    for j = 1:n
    like(j) = exp((lambda2(i)-lambda1(i))*j)*(lambda1(i)/lambda2(i))^sum(y(1:j));       
    end
    like = like/sum(like);
    theta(i) = randsample(1:n,1,true,like);   
end
b = 1;
tm = mean(theta(b:m));
lm2 = mean(lambda2(b:m));
lm1 = mean(lambda1(b:m));
subplot(3,1,1);
plot(lambda1);
ylabel('lambda1');
hold on;
subplot(3,1,2);
plot(lambda2);
ylabel('lambda2');
hold on;
subplot(3,1,3);
plot(theta);
ylabel('theta');
hold off;
figure();

subplot(3,1,1);
histogram(theta, 'Normalization','pdf');
ylabel('theta');
hold on;
ksdensity(theta);
hold on;
subplot(3,1,2);
histogram(lambda1, 'Normalization','pdf');
ylabel('lambda1');
hold on;
ksdensity(lambda1);
hold on;
subplot(3,1,3);
histogram(lambda2, 'Normalization','pdf');
ylabel('lambda2');
hold on;
ksdensity(lambda2);

subplot(3,1,1);
autocorr(lambda2,100,[],0)
title('for lambda2');
hold on;
subplot(3,1,2);
autocorr(lambda1,100,[],0)
title('for lambda1');
hold on;
subplot(3,1,3);
autocorr(theta,100,[],0)
title('for theta');



xi_L = prctile(lambda2, 2.5); 
xi_U = prctile(lambda2, 97.5); 
fprintf('(%f,%f)\n',xi_L, xi_U);

xi_L = prctile(lambda1, 2.5); 
xi_U = prctile(lambda1, 97.5); 
fprintf('(%f,%f)\n',xi_L, xi_U);

xi_L = prctile(theta, 2.5); 
xi_U = prctile(theta, 97.5); 
fprintf('(%f,%f)\n',xi_L, xi_U);


tu =unique(theta','rows');
l1u =unique(lambda1','rows');
l2u =unique(lambda2','rows');

min(theta);