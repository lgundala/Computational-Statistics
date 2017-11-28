clc;
clear;
lambda = 2;
n = 10000;
Z_i = zeros(1,n);
x_ = zeros(1,n);
x = zeros(1,n);
j=1;
ld = [2.2,4];
%%
%standard
for i = 1:n
    x_(i) = poissrnd(lambda);
Z = (x_(i) -2)/(sqrt(2/25));
if Z>=1.645
    x(j) = x_(i);
    Z_i(j) = Z;
    j= j+1;
end
end
est_Z = sum(Z_i)/j;
interval_m = [(1.645*sqrt(2/25) +2),(est_Z*sqrt(2/25) +2)];
%%
%antithetic
j = 1;
for i = 1:n/2
    x_(i) = poissrnd(lambda);
Z = (x_(i) -2)/(sqrt(2/25));
if Z>=1.645
    x(j) = x_(i);
    Z_i(j) = Z;
    j= j+1;
end
end
g = n/2;
for i = n/2+1:n
    x_(i) = -x_(g);
    g = g-1;
Z = (x_(i) -2)/(sqrt(2/25));
if Z>=1.645
    x(j) = x_(i);
    Z_i(j) = Z;
    j= j+1;
end
end
est_Z_a = sum(Z_i)/j;
interval_a = [(1.645*sqrt(2/25) +2),(est_Z_a*sqrt(2/25) +2)];
%%
%importance sampling
lambda_i = 2.4653;
x_i = zeros(1,n);
u = normrnd(0,1);
j= 1;
for i = 1:n
    x_(i) = poissrnd(lambda);
    x_i(i) = poissrnd(lambda_i);
y = x_(i)/(x_i(i)+10^-8);
Z = (x_(i) -2)/(sqrt(2/25));
if u<= y
    if Z>=1.645
    x(j) = x_(i);
    Z_i(j) =(x_(i) -2)/(sqrt(2/25));
    j= j+1;
    end
end
end
est_Z_i = sum(Z_i)/j;
interval_i = [(1.645*sqrt(2/25) +2),(est_Z_i*sqrt(2/25) +2)];

%% 
%importance sample with weights
lambda_i = 2.4653;
x_i = zeros(1,n);
u = normrnd(0,1);
j= 1;
w = zeros(1,n);
for i = 1:n
    x_(i) = poissrnd(lambda);
    x_i(i) = poissrnd(lambda_i);
y = x_i(i)/(x_(i));
Z = (x_(i) -2)/(sqrt(2/25));
if u<=y
    if Z>=1.645
    w(j) = y;
    x(j) = x_(i);
    Z_i(j) =(x_(i) -2)/(sqrt(2/25));
    j= j+1;
    end
end
end
w(isinf(w))=0;
j = j-1;
s = datasample(Z_i,n,'Weights',w);
est_Z_iw = sum(Z_i)/j;
interval_iw = [(1.645*sqrt(2/25) +2),(est_Z_iw*sqrt(2/25) +2)];
%%
% control variate
s_w = sum(w)/j;
Est_control_v = est_Z_iw -(s_w-1);
interval_cv = [(1.645*sqrt(2/25) +2),(Est_control_v*sqrt(2/25) +2)];

