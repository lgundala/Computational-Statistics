clc;
clear;

n = 10000;
m1 = 5000;
w1 = zeros(n,1);
s1= zeros(m1,1);
mw1 = zeros(m1,1);
mx1 = zeros(m1,1);
ew1 = zeros(n,1);
for i = 1: m1
x = rand(n,1);
q = exp((-abs(x).^3)./3);
g = 2*exp((-abs(x).^3));
w1 = g./q;
wi = w1./sum(w1);
for k = 1:n-1
ew1(k,:) = (x(k,:)*wi(k,:));
end
mw1(i,:) = sum(ew1);
mx1(i,:) = mean(x);
s1(i,:) = (std(ew1))^2;
end
sv = mean(s1);

%%
x1 = zeros(n,1);
likelihood = zeros(n,1);
e1 = zeros(n,1);
i=1;j=1;
while(j<=n)
 u=rand;
 x = normrnd(0,1);
 q = exp((-abs(x)^3)/3);
 g = 2*exp((-abs(x)^3));
 f = g/q;
 if (u<=f)
     x1(j,:) = x;
     likelihood(j,:) = q;
     w1(i,:) = g/q;
     e1(j,:) = g;
     j= j+1;
 end 
 i=i+1;
end
j=j-1;
wi = w1./sum(w1);
mw2 = sum(wi);
mx2 = mean(x);
% plot(x1,likelihood,'.');
% hold on;
% plot(x1,e1,'.');

%%
n1 = 10000;
m1 = 5000;
w1 = zeros(n1,1);
s1= zeros(m1,1);
mw1 = zeros(m1,1);
mx1 = zeros(m1,1);
ew1=zeros(n,1);
for i = 1: m1
x = rand(n1,1);
x1 = sort(x);
q = exp((-abs(x).^3)./3);
g = 2*exp((-abs(x).^3));
w1 = g./q;
wi = w1./sum(w1);
for k = 1:n-1
ew1(k,:) = (x(k+1,:)-x(k,:))*wi(k,:);
end
mw1(i,:) = sum(ew1);
mx1(i,:) = mean(x);
s1(i,:) = (std(ew1))^2;
end
em = mean(mw1);
ms = std(mw1);
s3 = mean(s1);


