clc;
clear;
%%
%SIR algorithm
A = importdata('coal.dat');
y = A.data(:,2);
n = length(y);
m = 5000;
a1 = 3;
a2 = 3;
l1 = zeros(1,1000);
l2 = zeros(1,1000);
lambda1 = zeros(1,m);
lambda2 = zeros(1,m);
theta = zeros(1,n);
like = zeros(1,n);
t = zeros(1,1000);
theta(1) = gamrnd(10,0.1,1,1);
lambda1(1) = 1;
lambda2(1) = 1;
b1 = 1;
b2 = 1;
c1 = 0;
c2 =0;
d1 = 1;
d2 = 1;
i =1;
for k = 2:m
     kk = theta(k-1);
    t = a1+sum(y(1:kk));
    lam = theta(k-1);
    lambda1(k) = gamrnd(t,1/lam,1,1);
    t = a2+sum(y) - sum(y(1:kk));
    lam = n-kk+b2;
    %lambda2(i) = gamrnd(t,1/lam,1,1);
    lambda2(k) = lambda1(k)*exp(unifrnd(log(1/8),log(2)));
    if((lambda1(k)/lambda2(k))>=unifrnd(0,1))
        l1(i) = lambda1(k);
        l2(i) = lambda2(k);
        t(i) = theta(k);
    for j = 1:n
    like(j) = exp((lambda2(k)-lambda1(k))*j)*(lambda1(k)/lambda2(k))^sum(y(1:j));       
    end
    like = like/sum(like);
    
    i = i+1;
    end
    b1 = gamrnd(a1+c1,1/(lambda1(i)+d1),1,1);
    b2 = gamrnd(a2+c2,1/(lambda1(i)+d2),1,1);
    theta(k) = randsample(1:n,1,true,like);
end
scatter(l1,l2)
b = 1;
tm = mean(t(b:m));
lm2 = mean(l2(b:m));
lm1 = mean(l1(b:m));
subplot(3,1,1);
plot(l1);
ylabel('lambda1');
hold on;
subplot(3,1,2);
plot(l2);
ylabel('lambda2');
hold on;
subplot(3,1,3);
plot(t);
ylabel('theta');
hold off;
figure();

subplot(3,1,1);
histogram(t, 'Normalization','pdf');
ylabel('theta');
hold on;
ksdensity(t);
hold on;
subplot(3,1,2);
histogram(l1, 'Normalization','pdf');
ylabel('lambda1');
hold on;
ksdensity(l1);
hold on;
subplot(3,1,3);
histogram(l2, 'Normalization','pdf');
ylabel('lambda2');
hold on;
ksdensity(l2);
figure();
subplot(3,1,1);
autocorr(l2,100,[],0)
title('for lambda2');
hold on;
subplot(3,1,2);
autocorr(l1,100,[],0)
title('for lambda1');
hold on;
subplot(3,1,3);
autocorr(theta,100,[],0)
title('for theta');
figure();

xi_L = prctile(l2, 2.5); 
xi_U = prctile(l2, 97.5); 
fprintf('(%f,%f)\n',xi_L, xi_U);
figure();
xi_L = prctile(l2, 2.5); 
xi_U = prctile(l2, 97.5); 
fprintf('(%f,%f)\n',xi_L, xi_U);
figure();
xi_L = prctile(t, 2.5); 
xi_U = prctile(t, 97.5); 
fprintf('(%f,%f)\n',xi_L, xi_U);


tu =unique(t','rows');
l1u =unique(l1','rows');
l2u =unique(l2','rows');