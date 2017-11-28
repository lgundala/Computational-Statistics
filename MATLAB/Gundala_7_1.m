clc;
clear;
%%
mu = [7,10];
delta = 0.7;
sigma = [0.5,0.5];
n = 200;
classlabel = binornd(1,delta,n,1);
y = classlabel.*normrnd(mu(1),sigma(1),[n,1])+(1-classlabel).*normrnd(mu(2),sigma(2),[n,1]);
histogram(y,20,'Normalization','pdf');
hold on;
ksdensity(y);
title('Histogram of 200 observations');
xlabel('y');
ylabel('Density');
%%
%Metropolis-Hastings Algorithm
likelihood = @(delta,y) prod(delta*normpdf(y,mu(1), sigma(1))+(1-delta)*normpdf(y,mu(2),sigma(2)));
nMC = 10000;
delta_i = rand;
deltas = zeros(nMC,1);
for i = 1:nMC
    delta_star = rand(1);
    MHratio = likelihood(delta_star,y)/likelihood(delta_i,y);
    if(rand<MHratio)
        delta_i = delta_star;
    end
    deltas(i) = delta_i;
end
% subplot(3,1,1);
plot(deltas);
xlabel('t');
ylabel('delta(t)');
title('MCMC');
% subplot(3,1,1);
% autocorr(deltas,40,[],0);
% title('MCMC');
%%
% random walk
nRW = nMC;
u = zeros(nRW,1);
u(1) = 2*rand(1,1)-1;
deltas = zeros(nRW,1);
likelihood = @(delta,y) prod(delta*normpdf(y,mu(1), sigma(1))+(1-delta)*normpdf(y,mu(2),sigma(2)));
delta_i = rand;
for i = 1:nRW
    delta_star = delta_i + u(i);
    RPratio = likelihood(delta_star,y)/likelihood(delta_i,y);
    if(rand<RPratio)
        delta_i = delta_star;
    end
    deltas(i) = delta_i;
    u(i+1) = 2*rand(1,1)-1;
    
end
% subplot(3,1,2);
plot(deltas);
xlabel('t');
ylabel('delta(t)');
title('Random Walk');
% subplot(3,1,2);
% autocorr(deltas,40,[],0);
% title('Random Walk');
%%
% parameterized random walk

nRP = nRW;
u = zeros(nRP,1);
u(1) = 2*rand(1,1)-1;
deltas = zeros(nRP,1);
deltas(1) = exp(u(1))/(1+exp(u(1)));
likelihood = @(delta,y) sum(log(delta*normpdf(y,mu(1), sigma(1))+(1-delta)*normpdf(y,mu(2),sigma(2))));
delta_i = rand;
for i = 1:nRP-1
    u(i+1) = u(i) + (2*rand(1,1)-1);
    deltas(i+1) = exp(u(i+1))/(1+exp(u(i+1)));
    RPratio = exp(likelihood(deltas(i+1),y) - likelihood(deltas(i),y))*exp(u(i))/exp(u(i+1));
    if(RPratio<1)
      if(binornd(1,RPratio)==0)
          deltas(i+1) = deltas(i);
          u(i+1) = u(i);
      end
    end 
end
% subplot(3,1,3);
plot(deltas);
xlabel('t');
ylabel('delta(t)');
title('Random Walk Re-Parameterized');

% subplot(3,1,3);
% autocorr(deltas,40,[],0);
% title('Random Walk Parameterized');
