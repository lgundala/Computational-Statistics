clear;
clc;
x = [8,3,4,3,1,7,2,6,2,7];
n=5000;
X = factorial(x);
x_1 = mean(x);
lambda1 = zeros(n,1);
j=1;
i=1;
lambda_old = zeros(n,1);
%% 
while(j<=n)
 u=rand;
 lambda= exp(normrnd(log(4),0.5));
 lambda_old(i,:) = lambda;
 likeQ = (lambda^(sum(x)))*exp(-10*lambda)/prod(X);
 likeE = (4.3^(sum(x)))*exp(-10*4.3)/prod(X);
 f = likeQ/likeE;
 if (u<=f)
     lambda1(j,:) = lambda;
     j= j+1;
 end 
 i=i+1;
end
j=j-1;

 L = mean(lambda1(1:j,:));
disp((j)*100/i);
la = (0.02:0.02:20)';
sx = sum(x);
px = prod(X);
lasx = la.^sx;
la10 = la*10;
a = lognpdf(la,log(4),0.5);
likeQ = a.*(lasx.*exp(-la10))/px;
likeE = a.*(4.3^sx)*exp(-10*4.3)/px;
hold on;
plot(la, likeQ,'.');
hold on;
plot(la, likeE);
xlabel('lambda');
ylabel('Unnmralized density');
legend('q(lambda|x)','e(lambda)');
g= figure();
hold off;
%%
[f,yout]=hist(lambda1,4);
bar(yout,f/sum(f));
hold on;
[t,u] = ksdensity(lambda1);
plot(u,t);
hold on;
fun = @(la) lognpdf(la,log(4),0.5).*((la.^sum(x)).*exp(-la*10))/prod(X);
funq = likeQ/integral(fun,0,Inf,'ArrayValued',true);
plot(la,funq);
legend('histogram','estimated positerior probability density curve','truepositerior probability density curve');