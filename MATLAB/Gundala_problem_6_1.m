clear;
clc;
r=2;
n=5000;
j=1;
x = zeros(n,1);
z = zeros(n,1);
a=r-(1/3);
b = 1/(sqrt(9*a));
u = zeros(n,1);
y = zeros(n,1);
i=1;
%% 
while(j<=n)
z_i = normrnd(0,1);
u(i,:)= rand;
disp(u(i,:));
out = a*(1+(b*z_i))^3;
y(i,:) = out;
if(out>0)
 z_1 = (z_i^2)/2;
out_1 = out/a;
l_o = a*log(out_1);
sum1 = (z_1 + l_o + a) - out;
f = exp(sum1);
disp(f);
if(u(i,:)<=f)
    z(j) = z_i;
    x(j) = out;
    j=j+1;
end
end
i=i+1;
end
j = j-1;
z1 = z.^2/2;
lo = a*log(x./a);
la = (0.02:0.02:20)';
ex = exp((lo+a)-x);
disp('acceptance rate');
disp(j*100/i);
plot(z,exp(-z1),'.');
hold on;
plot(z,ex,'.');
hold off;
legend('e','q');
f=figure;
%%
[e,xout]=hist(x,10);
bar(xout,e/sum(e));
hold on;
[o,p] = ksdensity(x);
plot(p,o);
hold on;
y = gampdf(la,r,1);
plot(la,y);
legend('histogram','estimated probability density curve','gamma(r,1)');