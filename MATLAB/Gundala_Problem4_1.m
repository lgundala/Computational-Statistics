clc;
clear;
nc = 85;
ni = 196;
l = 40;
pt_t1 = zeros(l,1);
pc_t1 = zeros(l,1);
pi_t1 = zeros(l,1);
pt_t = zeros(l,1);
pc_t = zeros(l,1);
pi_t = zeros(l,1);
r_t = zeros(l,1);
dc_t = zeros(l,1);
di_t = zeros(l,1);
dt_t = zeros(l,1);
nt = 341;
nu = 578;
n = ni+nc+nt;
pt_t(1,:) = 0.333333;
pc_t(1,:) = 0.333333;
pi_t(1,:) = 0.333333;
for i = 1 :l
ncc = nc*(pc_t(i,:))^2/((pc_t(i,:))^2 + 2*pc_t(i,:)*pi_t(i,:)+2*pc_t(i,:)*pt_t(i,:));
nci = 2*nc*pc_t(i,:)*pi_t(i,:)/((pc_t(i,:))^2 + 2*pc_t(i,:)*pi_t(i,:)+2*pc_t(i,:)*pt_t(i,:));
nct = 2*nc*pc_t(i,:)*pt_t(i,:)/((pc_t(i,:))^2 + 2*pc_t(i,:)*pi_t(i,:)+2*pc_t(i,:)*pt_t(i,:));
nii = ni*(pi_t(i,:))^2/(pi_t(i,:)^2 + 2*pi_t(i,:)*pt_t(i,:));
nit = 2*ni*pi_t(i,:)*pt_t(i,:)/(pi_t(i,:)^2 + 2*pi_t(i,:)*pt_t(i,:));
ntt = nt;
nii_u = nu*((pi_t(i,:))^2)/((pi_t(i,:)+pt_t(i,:))^2);
ntt_u = nu*((pt_t(i,:))^2)/((pi_t(i,:)+pt_t(i,:))^2);
nit_u = nu*(2*pi_t(i,:)*pt_t(i,:))/((pi_t(i,:)+pt_t(i,:))^2);
pc_t1(i,:) = ((2*ncc)+nci+nct)/(2*n);
pi_t1(i,:) = (2*nii+nit+nci+2*nii_u+nit_u)*(2*n - nci-nct)/((2*n-nci-nct+2*nu)*(2*n));
pt_t1(i,:) = (2*ntt+nit+nct+2*ntt_u +nit_u)*(2*n - nci-nct)/((2*n-nci-nct+2*nu)*(2*n));
pc_t(i+1,:) = pc_t1(i,:);
pt_t(i+1,:) = pt_t1(i,:);
pi_t(i+1,:) = pi_t1(i,:);
pi1 = pi_t(i,:);
pii = pi_t(i+1,:);
pt1 = pt_t(i,:);
ptt = pt_t(i+1,:);
if((abs(pi1-pii)<=10^-8)&&(abs(pt1-ptt)<=10^-8))
    break;
end

end
x = sprintf('values converged at %d iteration',i);
disp(x);
pc1 = pc_t(i,:);
pi1 = pi_t(i,:);
pt1 = pt_t(i,:);
p_t = [pc_t'; pi_t']';
l=i;

disp('     pc       pi       pt');
fprintf('%8.6f %8.6f %8.6f\n', [pc_t(1:l),pi_t(1:l),pt_t(1:l)]');
 