days = [0,8,28,41,63,79,97,117,135,154];%from table 2.3 time
beetles = [2,47,192,256,768,896,1120,896,1184,1024];%from table 2.3
l_beetles = log(beetles);
rt = ones(10,1);
c = ones(10,1);
theta = [rt c];
n_0 = l_beetles(1);
fn = @(days,r,k) ((k*n_0)./(n_0+(k-n_0)*exp(-r*days)));
g_theta = @(l_beetles,days,r,k) (1/2)*sum((l_beetles-fn).^2);
dlogfdr = @(days,r,k)((fn(days,r,k).^2).*days*(k-n_0).*exp(-r*days))./(k*n_0))./(fn(days,r,k));
dlogfdk = @(days,r,k)(((fn(days,r,k).^2).*(1-exp(-r*days))./(k^2)))./(fn(days,r,k));
dfdr = @(days,r,k)((fn(days,r,k).^2).*days*(k-n_0).*exp(-r*days))./(k*n_0);
dfdk = @(days,r,k)((fn(days,r,k).^2).*(1-exp(-r*days))./(k^2));
d2edr2 = @(days,r,k)(dfdr(days,r,k)).^2./(fn(days,r,k)).^2 +(((((fn(days,r,k).^2).*days.^2*(k-n_0).*exp(-r*days))./(k*n_0)).*(((2*fn(days,r,k)*(k-n_0).*exp(-r*days))/(k*n_0))-1))./fn(days,r,k));
d2edrk = @(days,r,k)(dfdr(days,r,k).*dfdk(days,r,k))./(fn(days,r,k)).^2 +((((fn(days,r,k).^2).*days.*exp(-r*days))./(k^2)).*((2*fn(days,r,k)*(k-n_0).*(1-exp(-r*days))/(k*n_0))+1))./(fn(days,r,k));
d2edk2 = @(days,r,k)(dfdk(days,r,k)).^2./(fn(days,r,k)).^2- ((-2*(fn(days,r,k).^2).*(1-exp(-r*days))./k^3).*((2*fn(days,r,k).*(1-exp(-r*days))./k)-1))./(fn(days,r,k)).^2;

for i = 1:10
    old_theta = theta(i,:);
    r = old_theta(1,1);
    k = old_theta(1,2);
    l_fn = log(fn(days, r, k));
    e = (l_beetles-l_fn)';
    jacob_theta = [(dlogfdr(days,r,k));(dlogfdk(days,r,k))]';
    g = -jacob_theta'*e;
    a11 = d2edr2(days,r,k);
    a12 = d2edrk(days,r,k);
    a21 = d2edrk(days,r,k);
    a22 = d2edk2(days,r,k);
    b11 = 0;
    b12 = 0;
    b21 = 0;
    b22 = 0;
    for k = 1:5
        b11 = b11+a11(k)*e(k);
        b12 = b12+a12(k)*e(k);
        b21 = b12;
        b22 = b22+a22(k)*e(k);
    end
    h = jacob_theta'*jacob_theta+[b11 b12;b21 b22];
    theta(i+1,:) = old_theta(1,:)-(((h)^-1)*g)';
end
newton = ((h)^-1);
e_theta1 = sqrt(newton(1,1));
e_theta2 = sqrt(newton(2,2));
correlation = newton(1,2)/(e_theta1*e_theta2);
disp(jacob_theta);
disp(h);
disp(newton);
disp(e_theta1);
disp(e_theta2);
disp(correlation);
