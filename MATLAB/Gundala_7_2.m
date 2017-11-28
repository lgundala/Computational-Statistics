clc;
clear;
%% initial values
%Metropolis-Hastings Algorithm
mu = [0,7,15];
delta = 0.7;
sigma = [0.01,0.01];
% sigma = [0.5,0.5];    % for new proposal distribution
n = 10000;
likelihood = @(x) prod(delta*(normrnd(mu(1),sigma(1)))+(1-delta)*(normrnd(mu(1),sigma(2))));
%%
pd = makedist('Normal','mu',mu(1),'sigma',sigma(1));
x_star = random(pd);
x = zeros(n,1);
x(1) = mu(1);
for i = 1:n-1
    R = (likelihood(x_star)*x(i))/(likelihood(x(i))*x_star);
    if R<=1
        x(i+1) = x_star;
    else
        x(i+1) = x(i);
    end
    pd = makedist('Normal','mu',x(i+1),'sigma',sigma(1));
    x_star = random(pd);
end
subplot(3,1,1);
% plot(x);
% autocorr(x,100,[],0);
histogram(x,20,'Normalization','pdf');
hold on;
ksdensity(x);
title('x(0)=0');

%%
pd = makedist('Normal','mu',mu(2),'sigma',sigma(1));
x_star = random(pd);
x = zeros(n,1);
x(1) = mu(2);
for i = 1:n-1
    R = (likelihood(x_star)*x(i))/(likelihood(x(i))*x_star);
    if R<=1
        x(i+1) = x_star;
    else
        x(i+1) = x(i);
    end
    pd = makedist('Normal','mu',x(i+1),'sigma',sigma(1));
    x_star = random(pd);
end
subplot(3,1,2);
% autocorr(x,100,[],0);
histogram(x,20,'Normalization','pdf');
hold on;
ksdensity(x);
title('x(0)=7');

%%
pd = makedist('Normal','mu',mu(3),'sigma',sigma(1));
x_star = random(pd);
x = zeros(n,1);
x(1) = mu(3);
for i = 1:n-1
    R = (likelihood(x_star)*x(i))/(likelihood(x(i))*x_star);
    if R<=1
        x(i+1) = x_star;
    else
        x(i+1) = x(i);
    end
    pd = makedist('Normal','mu',x(i+1),'sigma',sigma(1));
    x_star = random(pd);
end
subplot(3,1,3);
% autocorr(x,100,[],0);
histogram(x,20,'Normalization','pdf');
hold on;
ksdensity(x);
title('x(0)=15');

