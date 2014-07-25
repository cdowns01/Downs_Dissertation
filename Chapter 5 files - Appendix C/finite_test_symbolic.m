clear all
close all
t=cputime;
%Problem statement
% n''(x) = 12 x^3 over 0 to 1

%Set number of intervals to divide domain into
totalLength = 1;
numIntervals = 40;
intervalLength = totalLength/numIntervals;


%Define general basis functions for interval 0 to 1

%Phi function is chosen such that phi_right has f(0) = 1, f(1)=0, f'(0) = 1,
%f'(1)=0
%Phi function is chosen such that phi_left has f(0) = 0, f(1)=1, f'(0) = 1,
%f'(1)=0

%phi_right = 1-3.*x.^2+ 2.*x.^3;
%phi_left = 3*x.^2-2*x.^3; 
syms x

f1 =(1-3.*x.^2+ 2.*x.^3);
f2 =(3*x.^2-2*x.^3); 

%Psi function is chosen such that psi_right has f'(0) = 1, f'(1)=0, f(0) =
%0, f(1)=0
%Psi function is chosen such that psi_left has f'(0) = 0, f'(1)=1, f(0) = 0,
%f(1)=0
%psi_left = -(x.^2-x.^3);
%psi_right = x- 2*x.^2+x.^3;
f3 =-(x.^2-x.^3);
f4 =(x- 2*x.^2+x.^3);

%Create empty arrays to hold basis functions for each node
phi_right=sym(zeros(1,numIntervals+1));
phi_left=sym(zeros(1,numIntervals+1));
psi_right=sym(zeros(1,numIntervals+1));
psi_left=sym(zeros(1,numIntervals+1));

%Generate basis functions for each node (index-1)
for index=1:numIntervals+1
phi_right(index)= subs(f1,x,(x-(index-1)*intervalLength)/intervalLength);
phi_left(index)=subs(f2,x,(x-(index-2)*intervalLength)/intervalLength);

psi_right(index)=subs(f4,x,(numIntervals*x-(index-1)))*intervalLength;
psi_left(index)=subs(f3,x,(numIntervals*x-(index-2)))*intervalLength;
end

%Test to verify proper performance of basis functions at interval endpoints
% testPhR = zeros(1,numIntervals+1);
% testPhL = zeros(1,numIntervals+1);
% testPsR = zeros(1,numIntervals+1);
% testPsL = zeros(1,numIntervals+1);
% 
% for index = 1:numIntervals+1
%     testPhR(1,index) = subs(phi_right(index),x,(index-1)*intervalLength);
%     testPhR(2,index) = subs(phi_right(index),x,(index)*intervalLength);
%     testPhR(3,index) = subs(diff(phi_right(index),x),x,(index-1)*intervalLength);
%     testPhR(4,index) = subs(diff(phi_right(index),x),x,(index)*intervalLength);
%     
%     testPsR(1,index) = subs(psi_right(index),x,(index-1)*intervalLength);
%     testPsR(2,index) = subs(psi_right(index),x,(index)*intervalLength);
%     testPsR(3,index) = subs(diff(psi_right(index),x),x,(index-1)*intervalLength);
%     testPsR(4,index) = subs(diff(psi_right(index),x),x,(index)*intervalLength);
%         
%     testPhL(1,index) = subs(phi_left(index),x,(index-2)*intervalLength);
%     testPhL(2,index) = subs(phi_left(index),x,(index-1)*intervalLength);
%     testPhL(3,index) = subs(diff(phi_left(index),x),x,(index-2)*intervalLength);
%     testPhL(4,index) = subs(diff(phi_left(index),x),x,(index-1)*intervalLength);
%         
%     testPsL(1,index) = subs(psi_left(index),x,(index-2)*intervalLength);
%     testPsL(2,index) = subs(psi_left(index),x,(index-1)*intervalLength);
%     testPsL(3,index) = subs(diff(psi_left(index),x),x,(index-2)*intervalLength);
%     testPsL(4,index) = subs(diff(psi_left(index),x),x,(index-1)*intervalLength);
%     
%     
% end
%NOTE:At each node phi_left and psi_left are valid on the interval previous
%to the node, phi_right and psi_right are valid on the interval following
%the node.


x_axis = linspace(-1,1,1000);
testIndex=3;
plot(x_axis,subs(phi_right(testIndex),x,x_axis),x_axis,subs(phi_left(testIndex),x,x_axis),x_axis,subs(psi_right(testIndex),x,x_axis),x_axis,subs(psi_left(testIndex),x,x_axis))
axis([-0.2,1,-1,1])
%Legend('phi right','phi left', 'psi right', 'psi left')
%figure

%Set up empty matrix for holding coefficients
%Put Phi basis functions on odd rows, put Psi basis functions on even rows
%Subtract 2 rows and 2 columns to account for Boundary Condition
%limitations
matrix1 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);
% 
% %Create a list of coefficient variables for phi and psi functions
 Coeffs = sym(zeros(1, 2*(numIntervals+1)));
for k=1:numIntervals+1
     Coeffs(2*k-1) = sym(sprintf('c%d', k));
     Coeffs(2*k) = sym(sprintf('d%d', k));
end
 
%Fill center elements of array
%Let first term be summation basis function term, let second term be test
%function

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
     matrix1(2*k-4,2*k-2) = int(diff(phi_left(k),x,2)*phi_right(k-1),(k-2)*intervalLength, (k-1)*intervalLength);

    %Phi(k) Psi(k-1)
     matrix1(2*k-3,2*k-2) = int(diff(phi_left(k),x,2)*psi_right(k-1),(k-2)*intervalLength, (k-1)*intervalLength);
   
    %Phi(k) Phi(k)
    matrix1(2*k-2,2*k-2) = int(diff(phi_left(k),x,2)*phi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(diff(phi_right(k),x,2)*phi_right(k),(k-1)*intervalLength, (k)*intervalLength);
    
    %Phi(k) Psi(k)
    matrix1(2*k-1,2*k-2) = int(diff(phi_left(k),x,2)*psi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(diff(phi_right(k),x,2)*psi_right(k),(k-1)*intervalLength, (k)*intervalLength);
    
    %Phi(k) Phi(k+1)
    matrix1(2*k,2*k-2) = int(diff(phi_right(k),x,2)*phi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
    %Phi(k) Psi(k+1)
    matrix1(2*k+1,2*k-2) = int(diff(phi_right(k),x,2)*psi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix1(2*k-4,2*k-1) = int(diff(psi_left(k),x,2)*phi_right(k-1), (k-2)*intervalLength, (k-1)*intervalLength);

    %Psi(k) Psi(k-1)
    matrix1(2*k-3,2*k-1) = int(diff(psi_left(k),x,2)*psi_right(k-1), (k-2)*intervalLength, (k-1)*intervalLength);
   
    %Psi(k) Phi(k)
    matrix1(2*k-2,2*k-1) = int(diff(psi_left(k),x,2)*phi_left(k),(k-2)*intervalLength, (k-1)*intervalLength)+int(diff(psi_right(k),x,2)*phi_right(k),(k-1)*intervalLength, (k)*intervalLength);
    
    %Psi(k) Psi(k)
    matrix1(2*k-1,2*k-1) = int(diff(psi_left(k),x,2)*psi_left(k),(k-2)*intervalLength, (k-1)*intervalLength)+int(diff(psi_right(k),x,2)*psi_right(k),(k-1)*intervalLength,(k)*intervalLength);
    
    %Psi(k) Phi(k+1)
    matrix1(2*k,2*k-1) = int(diff(psi_right(k),x,2)*phi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
    %Psi(k) Psi(k+1)
    matrix1(2*k+1,2*k-1) = int(diff(psi_right(k),x,2)*psi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
end

%Handle 1st interval 'half case'
%Phis
%alpha is boundary conditions term, n'(0)=alpha*n(0)
alpha = 0;

matrix1(1,1) = int(diff(phi_right(1),x,2)*phi_right(1), 0, intervalLength)+alpha*int(diff(psi_right(1),x,2)*phi_right(1), 0, intervalLength);
%matrix1(2,1) = int(diff(phi_right(1),x,2)*psi_right(1), 0, intervalLength)+alpha*int(diff(psi_right(1),x,2)*psi_right(1), 0, intervalLength);
matrix1(2,1) = int(diff(phi_right(1),x,2)*phi_left(2),0, intervalLength)+alpha*int(diff(psi_right(1),x,2)*phi_left(2),0, intervalLength);
matrix1(3,1) = int(diff(phi_right(1),x,2)*psi_left(2),0, intervalLength)+alpha*int(diff(psi_right(1),x,2)*psi_left(2),0, intervalLength);


%Handle last interval 'half case'
%Phis
%beta is boundary condition term, n'(L) = beta n(L)
beta = 4/totalLength;

matrix1(2*(numIntervals)-2,2*(numIntervals+1)-2)= int(diff(phi_left(numIntervals+1),x,2)*phi_right(numIntervals),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength)+beta*int(diff(psi_left(numIntervals+1),x,2)*phi_right(numIntervals),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength);
matrix1(2*(numIntervals)-1,2*(numIntervals+1)-2)= int(diff(phi_left(numIntervals+1),x,2)*psi_right(numIntervals),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength)+beta*int(diff(psi_left(numIntervals+1),x,2)*psi_right(numIntervals),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength);
matrix1(2*(numIntervals),2*(numIntervals+1)-2)= int(diff(phi_left(numIntervals+1),x,2)*phi_left(numIntervals+1),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength)+beta*int(diff(psi_left(numIntervals+1),x,2)*phi_left(numIntervals+1),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength);
%matrix1(2*(numIntervals),2*(numIntervals+1)-2)= int(diff(phi_left(numIntervals+1),x,2)*psi_left(numIntervals+1),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength)+beta*int(diff(psi_left(numIntervals+1),x,2)*psi_left(numIntervals+1),(numIntervals-1)*intervalLength, (numIntervals)*intervalLength);    


%Handle interval 2 and n-1 case
%These intervals generate terms in the rows that would normally be in the
%for-loop matrix, but are omitted due to boundary conditions.  Due to the
%unusual number fo rows included in the final matrix, they must be input
%separately.

%Interval 2
%Phis
%Code is copied for code in for-loop above
k=2;
%Phi(k) Phi(k-1)
matrix1(2*k-3,2*k-2) = int(diff(phi_left(k),x,2)*phi_right(k-1),(k-2)*intervalLength, (k-1)*intervalLength);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
%matrix1(2*k-2,2*k-2) = int(diff(phi_left(k),x,2)*psi_right(k-1),(k-2)*intervalLength, (k-1)*intervalLength);
 
%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = int(diff(phi_left(k),x,2)*phi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(diff(phi_right(k),x,2)*phi_right(k),(k-1)*intervalLength, (k)*intervalLength);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = int(diff(phi_left(k),x,2)*psi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(diff(phi_right(k),x,2)*psi_right(k),(k-1)*intervalLength, (k)*intervalLength);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) = int(diff(phi_right(k),x,2)*phi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
%Phi(k) Psi(k+1)
matrix1(2*k+1,2*k-2) = int(diff(phi_right(k),x,2)*psi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-3,2*k-1) = int(diff(psi_left(k),x,2)*phi_right(k-1), (k-2)*intervalLength, (k-1)*intervalLength);

%Psi(k) Psi(k-1)
%matrix1(2*k-2,2*k-1) = int(diff(psi_left(k),x,2)*psi_right(k-1), (k-2)*intervalLength, (k-1)*intervalLength);
   
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = int(diff(psi_left(k),x,2)*phi_left(k),(k-2)*intervalLength, (k-1)*intervalLength)+int(diff(psi_right(k),x,2)*phi_right(k),(k-1)*intervalLength, (k)*intervalLength);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = int(diff(psi_left(k),x,2)*psi_left(k),(k-2)*intervalLength, (k-1)*intervalLength)+int(diff(psi_right(k),x,2)*psi_right(k),(k-1)*intervalLength,(k)*intervalLength);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) = int(diff(psi_right(k),x,2)*phi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
%Psi(k) Psi(k+1)
matrix1(2*k+1,2*k-1) = int(diff(psi_right(k),x,2)*psi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix1(2*k-4,2*k-2) = int(diff(phi_left(k),x,2)*phi_right(k-1),(k-2)*intervalLength, (k-1)*intervalLength);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix1(2*k-3,2*k-2) = int(diff(phi_left(k),x,2)*psi_right(k-1),(k-2)*intervalLength, (k-1)*intervalLength);
 
%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = int(diff(phi_left(k),x,2)*phi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(diff(phi_right(k),x,2)*phi_right(k),(k-1)*intervalLength, (k)*intervalLength);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = int(diff(phi_left(k),x,2)*psi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(diff(phi_right(k),x,2)*psi_right(k),(k-1)*intervalLength, (k)*intervalLength);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) = int(diff(phi_right(k),x,2)*phi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
%Phi(k) Psi(k+1)
%matrix1(2*k,2*k-2) = int(diff(phi_right(k),x,2)*psi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-4,2*k-1) = int(diff(psi_left(k),x,2)*phi_right(k-1), (k-2)*intervalLength, (k-1)*intervalLength);

%Psi(k) Psi(k-1)
matrix1(2*k-3,2*k-1) = int(diff(psi_left(k),x,2)*psi_right(k-1), (k-2)*intervalLength, (k-1)*intervalLength);
   
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = int(diff(psi_left(k),x,2)*phi_left(k),(k-2)*intervalLength, (k-1)*intervalLength)+int(diff(psi_right(k),x,2)*phi_right(k),(k-1)*intervalLength, (k)*intervalLength);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = int(diff(psi_left(k),x,2)*psi_left(k),(k-2)*intervalLength, (k-1)*intervalLength)+int(diff(psi_right(k),x,2)*psi_right(k),(k-1)*intervalLength,(k)*intervalLength);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) = int(diff(psi_right(k),x,2)*phi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);
    
%Psi(k) Psi(k+1)
%matrix1(2*k,2*k-1) = int(diff(psi_right(k),x,2)*psi_left(k+1),(k-1)*intervalLength, (k)*intervalLength);

%Multiply coefficients into LHS matrix
LHS=matrix1

RHS = 12*x^2;

%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrix(2*k-2) = int(RHS*phi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(RHS*phi_right(k),(k-1)*intervalLength,(k)*intervalLength);
    
    %Psi test function
    rhs_matrix(2*k-1) = int(RHS*psi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)+int(RHS*psi_right(k),(k-1)*intervalLength,(k)*intervalLength);
    
    
end
phi_left(k)
phi_right(k)
int(RHS*phi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)
int(RHS*phi_right(k),(k-1)*intervalLength,(k)*intervalLength)

int(RHS*psi_left(k),(k-2)*intervalLength,(k-1)*intervalLength)
int(RHS*psi_right(k),(k-1)*intervalLength,(k)*intervalLength)

%Handle 1st interval 'half case'
%Phi test function
rhs_matrix(1) = int(RHS*phi_right(1),0,intervalLength);
%Psi test function
%rhs_matrix(2) = int(RHS*psi_right(1),0,intervalLength);
    

%Handle last interval 'half case'
%Phi test function
rhs_matrix(2*(numIntervals+1)-2) = int(RHS*phi_left(numIntervals+1),(numIntervals-1)*intervalLength,(numIntervals)*intervalLength);

%Psi test function
%rhs_matrix(2*(numIntervals+1)) = int(RHS*psi_left(numIntervals+1),(numIntervals-1)*intervalLength,(numIntervals)*intervalLength);
    
RHS =rhs_matrix

%Solve simultaneous equations for basis function coefficients
Coeffs = LHS\RHS


guesstimate = zeros(2,numIntervals+1);
guesstimate(1,:)= linspace(0,1,numIntervals+1);
for index=2:length(guesstimate)
   guesstimate(2,index) = Coeffs(2*index-2);
end
guesstimate(2,1) = Coeffs(1);
guesstimate


y=linspace(0,1,100);

plot(y,y.^4)
hold on
plot(guesstimate(1,:),guesstimate(2,:))
hold off

time = cputime-t