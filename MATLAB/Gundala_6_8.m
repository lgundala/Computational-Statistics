clc;
clear;
n = 100000;
h = zeros(n,1);
t = zeros(1,n);
for i = 1:n
x = lognrnd(0,1);
e = normrnd(0,1);
tic;
y = exp(9 + 3*log(x) + e);
h(i,:) = y;
t(i) = toc;
end
ex = sum(h);
m = ex/n;
h_rao = zeros(n,1);
t_rao = zeros(1,n);
for i = 1:n
x_rao = lognrnd(0,1,[1,2]);
e = normrnd(0,1);
tic
y_rao = exp(9 + 3*log(x_rao(2)));
h_rao(i,:) = y_rao;
t_rao(i) = toc;
end
ex_rao = sum(h_rao);
m_rao = ex_rao/n;
plot(t,'*');
hold on;
plot(t_rao,'o');