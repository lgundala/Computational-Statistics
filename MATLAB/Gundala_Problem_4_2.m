clc;
clear;
nc = 85;
ni = 196;
n = zeros(8,1);
pt_t1 = zeros(8,1);
pc_t1 = zeros(8,1);
pi_t1 = zeros(8,1);
pt_t = zeros(8,1);
pc_t = zeros(8,1);
pi_t = zeros(8,1);
r_t = zeros(8,1);
dc_t = zeros(8,1);
di_t = zeros(8,1);
nt = 341;
pt_t(1,:) = 0.333333;
pc_t(1,:) = 0.333333;
pi_t(1,:) = 0.333333;
for i = 1 :9
ncc = nc*(pc_t(i,:))^2/((pc_t(i,:))^2 + 2*pc_t(i,:)*pi_t(i,:)+2*pc_t(i,:)*pt_t(i,:));
nci = 2*nc*pc_t(i,:)*pi_t(i,:)/((pc_t(i,:))^2 + 2*pc_t(i,:)*pi_t(i,:)+2*pc_t(i,:)*pt_t(i,:));
nct = 2*nc*pc_t(i,:)*pt_t(i,:)/((pc_t(i,:))^2 + 2*pc_t(i,:)*pi_t(i,:)+2*pc_t(i,:)*pt_t(i,:));
nii = ni*(pi_t(i,:))^2/(pi_t(i,:)^2 + 2*pi_t(i,:)*pt_t(i,:));
nit = 2*ni*pi_t(i,:)*pt_t(i,:)/(pi_t(i,:)^2 + 2*pi_t(i,:)*pt_t(i,:));
ntt = nt;
n(i,:) = nci+ncc+nii+ntt+nit+nct;
pc_t1(i,:) = ((2*ncc)+nci+nct)/(2*n(i,:));
pi_t1(i,:) = ((2*nii)+nci+nit)/(2*n(i,:));
pt_t1(i,:) = ((2*ntt)+nct+nit)/(2*n(i,:));
pc_t(i+1,:) = pc_t1(i,:);
pt_t(i+1,:) = pt_t1(i,:);
pi_t(i+1,:) = pi_t1(i,:);

end
pc1 = pc_t(i,:);
pi1 = pi_t(i,:);
pt1 = pt_t(i,:);
p_t = [pc_t'; pi_t']';

for i = 2:10
r_t(i,:) = norm(p_t(i,:) - p_t(i-1,:))/norm(p_t(i-1,:));
dc_t(i,:) = (pc_t(i,:)-pc1)/(pc_t(i-1,:) - pc1);
di_t(i,:) = (pi_t(i,:)-pi1)/(pi_t(i-1,:) - pi1);
end
disp('     pc       pi       R       Dc        Di');
fprintf('%8.6f %8.6f %8.6f %8.6f %8.6f\n', [pc_t(1:9),pi_t(1:9),r_t(1:9),dc_t(1:9), di_t(1:9)]');
 