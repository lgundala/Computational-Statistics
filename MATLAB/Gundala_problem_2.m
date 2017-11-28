gundala = importdata('baseball.dat');
salary = gundala.data(:,1);
y = gundala.data(:,2:28);
l_salary = log(salary);
mdl = stepwiselm(y,l_salary,linear);
