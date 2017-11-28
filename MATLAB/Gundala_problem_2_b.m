clc;
clear;
x = [162,267,271, 185, 111, 61, 27, 8, 3, 1];
lambda_t1 = zeros(10,1);
lambda_t2 = zeros(10,1);
pi_t = zeros(10,1);
fact = 1;
n = 10;
lambda_t1(1) = 1;
lambda_t2(1) = 1;
pi_t(1) = 1;
z_t = zeros(10,1);
for i= 1:n
    fact = fact*x(i);  
end
for k = 1: 100
pd1 = makedist('Poisson',lambda_t1(k,:));
y1 = pdf(pd1,x);
pd2 = makedist('Poisson',lambda_t2(k,:));
y2 = pdf(pd2,x);
f_x = (1-pi_t(k,:))*y1 + pi_t(k,:)*y2;
z_t = pi_t(k,:)*y2./(f_x+10^-8);
%disp(z_t);    
pi_t(k+1,:) = sum(z_t)/n;
lambda_t1(k+1,:) = sum((1-z_t)*x')/(n-sum(z_t));
lambda_t2(k+1,:) = sum(z_t*x')/sum(z_t);
t1 = pi_t(k,:);
t2 = pi_t(k+1,:);
if (abs(t1 - t2) <= 10^-8)
    disp(k);
    break;
end
end
y = [x',z_t']';
%disp(y);
disp('fitted probability');
disp(f_x);%fitted probability
disp('relative frequency');
disp(x'/n);
disp('    pi   lambda1    lambda2');
fprintf('%8.6f %8.6f %8.6f \n', [pi_t(1:k,:),lambda_t1(1:k,:),lambda_t2(1:k,:)]');