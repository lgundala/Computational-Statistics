clc;
clear;
%%
C = importdata('cancersurvival.dat');
ind1 = C.data(:,2)==1;
ind2 = C.data(:,2)==2;
S = C.data(ind1,1);
B = C.data(ind2,1);
s = log(S);
b = log(B);
N = 10000;
%%
% Stomach cancer
m = mean(s);
n = length(s);
m_star = zeros(N,1);
r_star = zeros(N,1);
std_Fhat = std(s);
for i = 1:N
    ind = randsample(n,n,'true');
    s_star = s(ind);
    m_star(i) = mean(s_star);
    std_star = sqrt(var(s_star));
    r_star(i) = (m_star(i)-m)/std_star;  
end

fprintf('for Stomach Cancer\n');

biashat = mean(m_star) - m; 
fprintf ('Bootstrap estimate of the bias is: %6.4f\n',  biashat)
fprintf ('Bias corrected bootstrap estimate is: %6.4f\n', m - biashat)
fprintf('standard error for corrected estimator is: %6.4f\n',std(m_star)/sqrt(length(m_star)))
xi_L = prctile(r_star, 2.5); 
xi_U = prctile(r_star, 97.5); 
fprintf ('95 percent bootstrap t C.I is: [ %f, %f ]\n', m - xi_U*std_Fhat, m - xi_L*std_Fhat)  
fprintf ('95 percent exp(C.I) using the simple bootstrap percentile method is: [ %f, %f ]\n', exp(prctile(m_star, 2.5)), exp(prctile(m_star, 97.5)))  
fprintf ('95 percent bootstrap t exp(C.I) is: [ %f, %f ]\n', exp(m - xi_U*std_Fhat), exp(m - xi_L*std_Fhat))

subplot(1,2,1)
histogram(r_star, 50); % Fig. 9.2
title('Frequency Plot of 10,000 vaules of R(xstar, Fhat)')
subplot(1,2,2)
histogram(r_star, 'Normalization','pdf'); % Fig. 9.2
title('Relative Frequency Plot')
%%
% Breast cancer
m = mean(b);
n = length(b);
m_star = zeros(N,1);
r_star = zeros(N,1);
std_Fhat = std(b);
for i = 1:N
    ind = randsample(n,n,'true');
    b_star = b(ind);
    m_star(i) = mean(b_star);
    std_star = sqrt(var(b_star));
    r_star(i) = (m_star(i)-m)/std_star;  
end


fprintf('\n');
fprintf('for Breast cancer\n');
biashat = mean(m_star) - m; 
fprintf ('Bootstrap estimate of the bias is: %6.4f\n',  biashat)
fprintf ('Bias corrected bootstrap estimate is: %6.4f\n', m - biashat)
fprintf('standard error for corrected estimator is: %6.4f\n',std(m_star)/sqrt(length(m_star)))
xi_L = prctile(r_star, 2.5); 
xi_U = prctile(r_star, 97.5); 
fprintf ('95 percent bootstrap t C.I is: [ %f, %f ]\n', m - xi_U*std_Fhat, m - xi_L*std_Fhat)  
fprintf ('95 percent exp(C.I) using the simple bootstrap percentile method is: [ %f, %f ]\n', exp(prctile(m_star, 2.5)), exp(prctile(m_star, 97.5)))
fprintf ('95 percent bootstrap t exp(C.I) is: [ %f, %f ]\n', exp(m - xi_U*std_Fhat), exp(m - xi_L*std_Fhat))
figure(1);
fprintf('for Stomach Cancer\n');
subplot(1,2,1)
histogram(r_star, 50); % Fig. 9.2
title('Frequency Plot of 10,000 vaules of R(xstar, Fhat)')
subplot(1,2,2)
histogram(r_star, 'Normalization','pdf'); % Fig. 9.2
title('Relative Frequency Plot')
%% for original data
% Stomach cancer
m = mean(S);
n = length(S);
m_star = zeros(N,1);
r_star = zeros(N,1);
std_Fhat = std(S);
for i = 1:N
    ind = randsample(n,n,'true');
    s_star = S(ind);
    m_star(i) = mean(s_star);
    std_star = sqrt(var(s_star));
    r_star(i) = (m_star(i)-m)/std_star;  
end
figure(2);
fprintf('for Stomach Cancer\n');
subplot(1,2,1)
histogram(r_star, 50); % Fig. 9.2
title('Frequency Plot of 10,000 vaules of R(xstar, Fhat)')
subplot(1,2,2)
histogram(r_star, 'Normalization','pdf'); % Fig. 9.2
title('Relative Frequency Plot')
biashat = mean(m_star) - m; 
fprintf ('Bootstrap estimate of the bias is: %6.4f\n',  biashat)
fprintf ('Bias corrected bootstrap estimate is: %6.4f\n', m - biashat)
fprintf('standard error for corrected estimator is: %6.4f\n',std(m_star)/sqrt(length(m_star)))
xi_L = prctile(r_star, 2.5); 
xi_U = prctile(r_star, 97.5); 
fprintf ('95 percent bootstrap t C.I is: [ %f, %f ]\n', m - xi_U*std_Fhat, m - xi_L*std_Fhat)  
fprintf ('95 percent C.I using the simple bootstrap percentile method is: [ %f, %f ]\n', prctile(m_star, 2.5), prctile(m_star, 97.5))
%%
% Breast cancer
m = mean(B);
n = length(B);
m_star = zeros(N,1);
r_star = zeros(N,1);
std_Fhat = std(B);
for i = 1:N
    ind = randsample(n,n,'true');
    b_star = B(ind);
    m_star(i) = mean(b_star);
    std_star = sqrt(var(b_star));
    r_star(i) = (m_star(i)-m)/std_star;  
end


fprintf('\n');
fprintf('for Breast cancer\n');
biashat = mean(m_star) - m; 
fprintf ('Bootstrap estimate of the bias is: %6.4f\n',  biashat)
fprintf ('Bias corrected bootstrap estimate is: %6.4f\n', m - biashat)
fprintf('standard error for corrected estimator is: %6.4f\n',std(m_star)/sqrt(length(m_star)))
xi_L = prctile(r_star, 2.5); 
xi_U = prctile(r_star, 97.5); 
fprintf ('95 percent bootstrap t C.I is: [ %f, %f ]\n', m - xi_U*std_Fhat, m - xi_L*std_Fhat)  
fprintf ('95 percent C.I using the simple bootstrap percentile method is: [ %f, %f ]\n', prctile(m_star, 2.5), prctile(m_star, 97.5))

figure(2);
fprintf('for Stomach Cancer\n');
subplot(1,2,1)
histogram(r_star, 50); % Fig. 9.2
title('Frequency Plot of 10,000 vaules of R(xstar, Fhat)')
subplot(1,2,2)
histogram(r_star, 'Normalization','pdf'); % Fig. 9.2
title('Relative Frequency Plot')