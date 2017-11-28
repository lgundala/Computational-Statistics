clc;
clear;
n = 100000;
s = zeros(1,n);
s_ = zeros(1,n);
s_0 = 50;
C = zeros(1,n);
K = 52;
r =0.05;
A = zeros(1,n);
sigma_0 = 0.5;
N = 30;
s_bar = zeros(1,n);
s_bar_1 = zeros(1,n);

x = normrnd(0,1);
p = normcdf([-1,1]);
z = p(2)-p(1);
T = 30;
Antithetic_A = zeros(1,n);
z_i = zeros(1,30);
theta = zeros(1,30);
theta_m = zeros(1,n);
s(1) = s_0*exp((r-(sigma_0^2))/365 + (sigma_0*z/sqrt(365)));
s_(1) = s_0*exp((r-(sigma_0^2))/365 + (sigma_0*z/sqrt(365)));
for i = 1:n
     e = 30;
for k = 1:e/2
    z_i(k) = normrnd(0,1);
end
g = e/2;
for k = e/2+1:e
    z_i(k) = -z_i(g);
    g = g-1;
end
for j = 1:29
    z = normrnd(0,1);
   s(j+1) = s(j) * exp((r-(sigma_0^2))/365 + (sigma_0*z/sqrt(365)));
   s_(j+1) = s_(j) * exp((r-(sigma_0^2))/365 + (sigma_0*z_i(j)/sqrt(365)));
  
end
for v = 1:30
c_3 = 1+ 1/N;
c_2 = sigma_0*sqrt((c_3*v/1095)*(1+(1/2*N)));
c_1 = (1/c_2)*(log(s_0/K)+((c_3*v/730)*(r - ((sigma_0^2)/2)))+((c_3*(sigma_0^2)*v/1095)*(1+(1/2*N))));
phi = normcdf(c_1);
phi_2 = normcdf(c_1-c_2);
theta(v) = s_0*phi*exp(-v*(r+(c_3*(sigma_0^2)/6))*((1-1/N)/730)) - (K*phi_2*exp(-r*v/365));
end
theta_m(i) = sum(theta)/T;
s_bar(i) = sum(s)/T;
s_bar_1(i) = sum(s_)/T;
  Antithetic_A(i) = exp(-(r*T)/365)* max(0,s_bar_1(i)-K);
     z = normrnd(0,1);
     s_T = s_0*exp((r-((sigma_0^2)/2))*T/365 + (sigma_0*z*sqrt(T/365)));
   C(i) = exp(-(r*T)/365)*max(0,s_T-K);
   A(i) = exp(-(r*T)/365)* max(0,s_bar(i)-K);
end
fair_price = sum(C)/n;
Asian_price = sum(A)/n;
anthetic_fair_price = sum(Antithetic_A)/n;
control_variate_theta = sum(theta_m)/n;
plot(A'.');
hold on;
plot(Antithetic_A'.');
hold on;
plot(theta_m,'.');


