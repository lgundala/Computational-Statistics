clc;
clear;
%%
% 4.2.b
N = 1500;
% al = [0.2367108,0.18435487,0.17445448,0.17214260,0.171290563,0.170665118,0.170094226,0.169588178,0.16916086,0.168810776,0.168528271];
n0 = 379;
n = [379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1]';
alpha = zeros(17,1);
beta = zeros(17,1);
mu = zeros(17,1);
lambda = zeros(17,1);
alpha(1) = 0.1;
beta(1) = 0.1;
mu(1) = 1;
lambda(1) = 1;
phi = zeros(17,1);
a = 0.1;
b = 0.1;
m = 1;
l = 1;
phi0 = @(a,b,m,l) a + b*poisspdf(0,m) + (1-a-b)*poisspdf(0,l);
z0 = a/phi0(a,b,m,l);
phtheta = @(a,b,m,l,i) b*poisspdf(i,m)*factorial(i) + (1-a-b)*poisspdf(i,l)*factorial(i);
ttheta = @(b,m,pha,i) b*poisspdf(i,m)*factorial(i)/pha;
ztheta = @(a,ph0) a/ph0;
ptheta= @(a,b,l,pha,i) (1-a-b)*poisspdf(i,l)*factorial(i)/pha;
t = zeros(17,1);
p = zeros(17,1);
p(1) = 1;
z=z0;
for j = 1:100
    
for k =1:17
phi(k) = phtheta(alpha(j),beta(j),mu(j),lambda(j),k-1);
t(k) = ttheta(beta(j),mu(j),phi(k),k-1);
z = ztheta(alpha(j),phi0(alpha(j),beta(j),mu(j),lambda(j)));
p(k) = ptheta(alpha(j),beta(j),lambda(j),phi(k),k-1);
end
t0 = ttheta(beta(1),mu(1),phtheta(alpha(1),beta(1),mu(1),lambda(1),0),0);
p0 = ptheta(alpha(1),beta(1),lambda(1),phtheta(alpha(1),beta(1),mu(1),lambda(1),0),0);
alpha(j+1) = n0*z/N;
beta(j+1) = (sum(n.*t))/N;
mu1 = 0;
l1 = 0;
for o=1:16
    mu1 = mu1 +o*n(o)*t(o);
    l1 = l1 +o*n(o)*p(o);
end
mu(j+1) = mu1/(sum(n.*t)+n0*t0);
lambda(j+1) = l1/(sum(n.*p)+n0*p0);
if(abs(alpha(j)-alpha(j+1))<10^-6)
    break;
end
end
fprintf('alpha\t\tbeta\t\tmu\t\tlambda\tR(relative convergence)\n');
for u = j-10:j+1
    fprintf('%f\t%f\t%f\t%f\t%f\n',alpha(u),beta(u),mu(u),lambda(u),(abs((alpha(u)-alpha(u-1))/alpha(u-1))));
end