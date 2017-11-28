% MATLAB code for the examples on Ch9 of
% Computational Statistics, 2nd edition, Givens & Hoeting,

x = [0.01 0.48 0.71 0.95 1.19 0.01 0.48 1.44 0.71 1.96 0.01 1.44 1.96];
y = [127.6 124.0 110.8 103.9 101.5 130.1 122.0 92.3 113.1 83.7 128.0 91.4 86.2];

n = length(x);
B = 10000;

beta_covb = regstats(y, x,'linear',{'beta', 'covb'});
%p = polyfit(x,y,1);
%thetahat = p(1)/p(2);
beta = beta_covb.beta;
covb = beta_covb.covb;

thetahat = beta(2)/beta(1);

std_Fhat = abs(thetahat)*sqrt(covb(2,2)/(beta(2)^2) + covb(1,1)/(beta(1))^2 - 2*covb(1,2)/(beta(1)*beta(2))); 
% see equation(9.16) on page 297

thetahat_star = zeros(B,1);
R_star = zeros(B,1);

for b = 1:B
    ind = randsample(n, n, 'true')';
    y_star = y(ind);
    x_star = x(ind);
    %p_star = polyfit(x_star,y_star,1);
    %thetahat_star(b) = p_star(1)/p_star(2);
    beta_covb_star = regstats(y_star,x_star,'linear',{'beta', 'covb'});
    beta_star = beta_covb_star.beta;
    covb_star = beta_covb_star.covb;
    thetahat_star(b) = beta_star(2)/beta_star(1);
    std_Fhat_star = sqrt(thetahat_star(b)^2*(covb_star(2,2)/beta_star(2)^2 + covb_star(1,1)/beta_star(1)^2 - 2*covb_star(1,2)/(beta_star(1)*beta_star(2))));
    R_star(b) = (thetahat_star(b) - thetahat)/std_Fhat_star;
    fprintf(' b = %g \n', b) 
end

% Example 9.3
figure(1)
subplot(1,2,1)
histogram(thetahat_star, 50); % Fig. 9.1
title('Fig 9.1(a) - Frequency Plot of 10,000 bootstrap estimates of beta_1/beta_0')
subplot(1,2,2)
histogram(thetahat_star,'Normalization','pdf')
title('Fig 9.1(b) - Relative Frequency Plot')

% Example 9.4
biashat = mean(thetahat_star) - thetahat; 
fprintf ('Bootstrap estimate of the bias is: %6.4f\n',  biashat)
fprintf ('Bias corrected bootstrap estimate is: %6.4f\n', thetahat - biashat)

% Example 9.5
fprintf ('95 percent C.I using the simple bootstrap percentile method is: [ %f, %f ]\n', prctile(thetahat_star, 2.5), prctile(thetahat_star, 97.5))  

% Example 9.7
xi_L = prctile(R_star, 2.5); 
xi_U = prctile(R_star, 97.5); 
fprintf ('95 percent bootstrap t C.I is: [ %f, %f ]\n', thetahat + xi_L*std_Fhat, thetahat + xi_U*std_Fhat)  
figure(2)
subplot(1,2,1)
histogram(R_star, 50); % Fig. 9.2
title('Fig 9.2(a) - Frequency Plot of 10,000 vaules of R(xstar, Fhat)')
subplot(1,2,2)
histogram(R_star, 'Normalization','pdf'); % Fig. 9.2
title('Fig 9.2(b) - Relative Frequency Plot')

