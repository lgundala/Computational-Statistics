clc;
clear;
%%
% 9.4 a
A = importdata('salmon.dat');
r = A.data(:,2);
s = A.data(:,3);
n = length(r);
B = 10000;
R = r.^-1;
S = s.^-1;
% fitvars = polyfit(S, R, 1);
% beta1 = fitvars(2);
% beta2 = fitvars(1);
beta_covb = regstats(R, S,'linear',{'beta', 'covb'});
beta = beta_covb.beta;
covb = beta_covb.covb;

thetahat = beta(1)/(1-beta(2));
std_Fhat = abs(thetahat)*sqrt(covb(2,2)/(beta(2)^2) + covb(1,1)/(beta(1))^2 - 2*covb(1,2)/(beta(1)*beta(2)));

thetahat_star = zeros(B,1);
r_star = zeros(B,1);

plot(r,s,'.');
for b = 1:B
    ind = randsample(n, n, 'true')';
    R_star = R(ind);
    S_star = S(ind);
    %p_star = polyfit(x_star,y_star,1);
    %thetahat_star(b) = p_star(1)/p_star(2);
    beta_covb_star = regstats(R_star,S_star,'linear',{'beta', 'covb'});
    beta_star = beta_covb_star.beta;
    covb_star = beta_covb_star.covb;
    thetahat_star(b) = beta_star(1)/(1-beta_star(2));
    std_Fhat_star = sqrt(thetahat_star(b)^2*(covb_star(2,2)/beta_star(2)^2 + covb_star(1,1)/beta_star(1)^2 - 2*covb_star(1,2)/(beta_star(1)*beta_star(2))));
    r_star(b) = (thetahat_star(b) - thetahat)/std_Fhat_star;
%     fprintf(' b = %g \n', b) 
end
fprintf('for cases\n');
figure(1)
subplot(1,2,1)
histogram(thetahat_star, 50);
title('Frequency Plot of 10,000 bootstrap estimates')
subplot(1,2,2)
histogram(thetahat_star,'Normalization','pdf')
title(' Relative Frequency Plot')

% Example 9.4
biashat = mean(thetahat_star) - thetahat; 
fprintf ('Bootstrap estimate of the bias is: %6.8f\n',  biashat)
fprintf ('Bias corrected bootstrap estimate is: %6.8f\n', thetahat - biashat)
fprintf('standard error for corrected estimator is: %6.8f\n',std(thetahat_star)/sqrt(length(thetahat_star)))
% Example 9.5
fprintf ('95 percent C.I using the simple bootstrap percentile method is: [ %f, %f ]\n', prctile(thetahat_star, 2.5), prctile(thetahat_star, 97.5))  

% Example 9.7
xi_L = prctile(r_star, 2.5); 
xi_U = prctile(r_star, 97.5); 
fprintf ('95 percent bootstrap t C.I is: [ %f, %f ]\n', thetahat - xi_U*std_Fhat, thetahat - xi_L*std_Fhat)  
figure(2)
subplot(1,2,1)
histogram(r_star, 50); % Fig. 9.2
title('Frequency Plot of 10,000 vaules of R(xstar, Fhat)')
subplot(1,2,2)
histogram(r_star, 'Normalization','pdf'); % Fig. 9.2
title('Relative Frequency Plot')
fprintf('  \n  ');
%%
%residuals
p = polyfit(R,S,1);
yfit = polyval(p,S);
yresid = R - yfit;
beta_covb = regstats(yfit, S,'linear',{'beta', 'covb'});
beta = beta_covb.beta;
covb = beta_covb.covb;

thetahat = beta(1)/(1-beta(2));
std_Fhat = abs(thetahat)*sqrt(covb(2,2)/(beta(2)^2) + covb(1,1)/(beta(1))^2 - 2*covb(1,2)/(beta(1)*beta(2)));


for b = 1:B
    ind = randsample(n, n, 'true')';
    y_bar = yresid(ind);
    y_true = yfit(ind) + yresid(ind);
    S_star = S(ind);
    %p_star = polyfit(x_star,y_star,1);
    %thetahat_star(b) = p_star(1)/p_star(2);
    beta_covb_star = regstats(y_true,S_star,'linear',{'beta', 'covb'});
    beta_star = beta_covb_star.beta;
    covb_star = beta_covb_star.covb;
    thetahat_star(b) = beta_star(1)/(1-beta_star(2));
    std_Fhat_star = sqrt(thetahat_star(b)^2*(covb_star(2,2)/beta_star(2)^2 + covb_star(1,1)/beta_star(1)^2 - 2*covb_star(1,2)/(beta_star(1)*beta_star(2))));
    r_star(b) = (thetahat_star(b) - thetahat)/std_Fhat_star;
%     fprintf(' b = %g \n', b) 
end
fprintf('     \n');
fprintf('for residuals\n');
figure(3)
subplot(1,2,1)
histogram(thetahat_star, 50); % Fig. 9.1
title('Frequency Plot of 10,000 bootstrap estimates')
subplot(1,2,2)
histogram(thetahat_star,'Normalization','pdf')
title('Relative Frequency Plot')
% Example 9.4
biashat = mean(thetahat_star) - thetahat; 
fprintf ('Bootstrap estimate of the bias is: %6.8f\n',  biashat)
fprintf ('Bias corrected bootstrap estimate is: %6.8f\n', thetahat - biashat)
fprintf('standard error for corrected estimator is: %6.8f\n',std(thetahat_star)/sqrt(length(thetahat_star)))
% Example 9.5
fprintf ('95 percent C.I using the simple bootstrap percentile method is: [ %f, %f ]\n', prctile(thetahat_star, 2.5), prctile(thetahat_star, 97.5))  

% Example 9.7
xi_L = prctile(r_star, 2.5); 
xi_U = prctile(r_star, 97.5); 
fprintf ('95 percent bootstrap t C.I is: [ %f, %f ]\n', (thetahat - xi_U*std_Fhat), (thetahat - xi_L*std_Fhat))  
figure(4)
subplot(1,2,1)
histogram(r_star, 50); % Fig. 9.2
title('Frequency Plot of 10,000 vaules of R(xstar, Fhat)')
subplot(1,2,2)
histogram(r_star, 'Normalization','pdf'); % Fig. 9.2
title('Relative Frequency Plot')
