clear all
close all
t=cputime;

%Version 1.1
%Calculates the approximate solution of a known test function
%Uses numerical integration instead of symbolic integration
%The length of the testing regions can be varied
%The number of intervals can be varied
%The starting point of the testing region can be varied

%Input the start and end points of the testing region and determine its
%length
startPoint = 0;
endPoint = 1;
totalLength = endPoint-startPoint;

%Set number of intervals to divide domain into
numIntervals = 10;
intervalLength = totalLength/numIntervals;

syms x

%RHS = n''(x)
%Corresponds to n(x) = x^4
RHS = 12*x^2;

%alpha is boundary conditions term, n'(0)=alpha*n(0)
if startPoint==0
    alpha = 0;
else
    alpha = 4/startPoint;
end

%beta is boundary condition term, n'(L) = beta n(L)
if endPoint==0
    beta = 0;
else
    beta = 4/endPoint;
end

%Define general basis functions for interval 0 to 1

%Phi function is chosen such that phi_right (f1) has f(0) = 1, f(1)=0, f'(0) = 1,
%f'(1)=0
%Phi function is chosen such that phi_left (f2) has f(0) = 0, f(1)=1, f'(0) = 1,
%f'(1)=0

f1 =(1-3.*x.^2+ 2.*x.^3);
f2 =(3*x.^2-2*x.^3); 

%Psi function is chosen such that psi_right (f4) has f'(0) = 1, f'(1)=0, f(0) =
%0, f(1)=0
%Psi function is chosen such that psi_left (f3) has f'(0) = 0, f'(1)=1, f(0) = 0,
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
phi_right(index)= subs(f1,x,(x-(index-1)*intervalLength-startPoint)/intervalLength);
phi_left(index)=subs(f2,x,(x-(index-2)*intervalLength-startPoint)/intervalLength);

psi_right(index)=subs(f4,x,(x-(index-1)*intervalLength-startPoint)/intervalLength)*intervalLength;%(numIntervals*x-(index-1)))*intervalLength;
psi_left(index)=subs(f3,x,(x-(index-2)*intervalLength-startPoint)/intervalLength)*intervalLength;%(numIntervals*x-(index-2)))*intervalLength;
end
%NOTE:At each node phi_left and psi_left are valid on the interval previous
%to the node, phi_right and psi_right are valid on the interval following
%the node.

%Test to verify proper performance of basis functions at interval endpoints
% testPhR = zeros(1,numIntervals+1);
% testPhL = zeros(1,numIntervals+1);
% testPsR = zeros(1,numIntervals+1);
% testPsL = zeros(1,numIntervals+1);
% 
% for index = 1:numIntervals+1
%     testPhR(1,index) = subs(phi_right(index),x,(index-1)*intervalLength+startPoint);
%     testPhR(2,index) = subs(phi_right(index),x,(index)*intervalLength+startPoint);
%     testPhR(3,index) = subs(diff(phi_right(index),x),x,(index-1)*intervalLength+startPoint);
%     testPhR(4,index) = subs(diff(phi_right(index),x),x,(index)*intervalLength+startPoint);
%     
%     testPsR(1,index) = subs(psi_right(index),x,(index-1)*intervalLength+startPoint);
%     testPsR(2,index) = subs(psi_right(index),x,(index)*intervalLength+startPoint);
%     testPsR(3,index) = subs(diff(psi_right(index),x),x,(index-1)*intervalLength+startPoint);
%     testPsR(4,index) = subs(diff(psi_right(index),x),x,(index)*intervalLength+startPoint);
%         
%     testPhL(1,index) = subs(phi_left(index),x,(index-2)*intervalLength+startPoint);
%     testPhL(2,index) = subs(phi_left(index),x,(index-1)*intervalLength+startPoint);
%     testPhL(3,index) = subs(diff(phi_left(index),x),x,(index-2)*intervalLength+startPoint);
%     testPhL(4,index) = subs(diff(phi_left(index),x),x,(index-1)*intervalLength+startPoint);
%         
%     testPsL(1,index) =
%     subs(psi_left(index),x,(index-2)*intervalLength+startPoint);
%     testPsL(2,index) = subs(psi_left(index),x,(index-1)*intervalLength+startPoint);
%     testPsL(3,index) = subs(diff(psi_left(index),x),x,(index-2)*intervalLength+startPoint);
%     testPsL(4,index) = subs(diff(psi_left(index),x),x,(index-1)*intervalLength+startPoint);
%     
%     
% end

 x_axis = linspace(-totalLength,totalLength,1000);
 testIndex=1;
 plot(x_axis,subs(phi_right(testIndex),x,x_axis),x_axis,subs(phi_left(testIndex),x,x_axis),x_axis,subs(psi_right(testIndex),x,x_axis),x_axis,subs(psi_left(testIndex),x,x_axis))
 axis([-totalLength,totalLength,-1,1])
 figure

%Set up empty matrix for holding coefficients
%Put Phi basis functions on odd rows, put Psi basis functions on even rows
%Subtract 2 rows and 2 columns to account for Boundary Condition
%limitations
matrix1 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);
 
%Create a list of coefficient variables for phi and psi functions

%Fill center elements of array
%Let first term be summation basis function term, let second term be test
%function

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix1(2*k-4,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Phi(k) Psi(k-1)
    matrix1(2*k-3,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Phi(k) Phi(k)
    matrix1(2*k-2,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(phi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k)
    matrix1(2*k-1,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(phi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Phi(k+1)
    matrix1(2*k,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k+1)
    matrix1(2*k+1,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix1(2*k-4,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Psi(k) Psi(k-1)
    matrix1(2*k-3,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Psi(k) Phi(k)
    matrix1(2*k-2,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k)
    matrix1(2*k-1,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi(k) Phi(k+1)
    matrix1(2*k,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k+1)
    matrix1(2*k+1,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
end

%Handle 1st node 'half case'
%Phis
matrix1(1,1) = quad(matlabFunction(diff(phi_right(1),x,2)*phi_right(1)), startPoint, intervalLength+startPoint)+alpha*quad(matlabFunction(diff(psi_right(1),x,2)*phi_right(1)),startPoint, intervalLength+startPoint);
matrix1(2,1) = quad(matlabFunction(diff(phi_right(1),x,2)*phi_left(2)),startPoint, intervalLength+startPoint)+alpha*quad(matlabFunction(diff(psi_right(1),x,2)*phi_left(2)),startPoint, intervalLength+startPoint);
matrix1(3,1) = quad(matlabFunction(diff(phi_right(1),x,2)*psi_left(2)),startPoint, intervalLength+startPoint)+alpha*quad(matlabFunction(diff(psi_right(1),x,2)*psi_left(2)),startPoint, intervalLength+startPoint);

%Handle last node 'half case'
%Phis
matrix1(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(matlabFunction(diff(phi_left(numIntervals+1),x,2)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(matlabFunction(diff(psi_left(numIntervals+1),x,2)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrix1(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(matlabFunction(diff(phi_left(numIntervals+1),x,2)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(matlabFunction(diff(psi_left(numIntervals+1),x,2)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrix1(2*(numIntervals),2*(numIntervals+1)-2)= quad(matlabFunction(diff(phi_left(numIntervals+1),x,2)*phi_left(numIntervals+1)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(matlabFunction(diff(psi_left(numIntervals+1),x,2)*phi_left(numIntervals+1)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);

%Handle node 2 and n-1 case
%These nodes generate terms in the rows that would normally be in the
%for-loop matrix, but are omitted due to boundary conditions.  Due to the
%unusual number fo rows included in the final matrix, they must be input
%separately.

%Node 2
%Phis
%Code is copied for code in for-loop above
k=2;
%Phi(k) Phi(k-1)
matrix1(2*k-3,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(phi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(phi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k+1)
matrix1(2*k+1,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-3,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
  
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Psi(k) Psi(k+1)
matrix1(2*k+1,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix1(2*k-4,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix1(2*k-3,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(phi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(phi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-4,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

%Psi(k) Psi(k-1)
matrix1(2*k-3,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    

%Multiply coefficients into LHS matrix
LHS=matrix1

%Calculates for RHS of equation (forcing function time test functions)

%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrix(2*k-2) = quad(matlabFunction(RHS*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(RHS*phi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi test function
    rhs_matrix(2*k-1) = quad(matlabFunction(RHS*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(RHS*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    
end

%Handle 1st interval 'half case'
%Phi test function
rhs_matrix(1) = quad(matlabFunction(RHS*phi_right(1)),startPoint,intervalLength+startPoint);    

%Handle last interval 'half case'
%Phi test function
rhs_matrix(2*(numIntervals+1)-2) = quad(matlabFunction(RHS*phi_left(numIntervals+1)),(numIntervals-1)*intervalLength+startPoint,(numIntervals)*intervalLength+startPoint);
   
RHSOut =rhs_matrix

%Solve simultaneous equations for basis function coefficients
Coeffs = LHS\RHSOut


%Collect Cj's in an array for display purposes
guesstimate = zeros(2,numIntervals+1);
guesstimate(1,:)= linspace(startPoint,totalLength+startPoint,numIntervals+1);
for index=2:length(guesstimate)
   guesstimate(2,index) = Coeffs(2*index-2);
end
guesstimate(2,1) = Coeffs(1)


%Plotting
y=linspace(startPoint,totalLength+startPoint,100);

plot(y,y.^4)
hold on
plot(guesstimate(1,:),guesstimate(2,:))
hold off

%Displays runtime of m-file
time = cputime-t