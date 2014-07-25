%Problem statement
% n''(x) = 12 x^3 over 0 to 1

%Set number of intervals to divide domain into
totalLength = 1;
numIntervals = 5;
intervalLength = totalLength/numIntervals;


%Define general basis functions for interval 0 to 1

%Phi function is chosen such that phi_right has f(0) = 1, f(1)=0, f'(0) = 1,
%f'(1)=0
%Phi function is chosen such that phi_left has f(0) = 0, f(1)=1, f'(0) = 1,
%f'(1)=0
syms x
%phi_right = 1-3.*x.^2+ 2.*x.^3;
%phi_left = 3*x.^2-2*x.^3; 

f1 = 1-3.*x.^2+ 2.*x.^3;
f2 = 3*x.^2-2*x.^3; 

%Psi function is chosen such that psi_right has f'(0) = 1, f'(1)=0, f(0) =
%0, f(1)=0
%Psi function is chosen such that psi_left has f'(0) = 0, f'(1)=1, f(0) = 0,
%f(1)=0
%psi_left = -(x.^2-x.^3);
%psi_right = x- 2*x.^2+x.^3;
f3 = -(x.^2-x.^3);
f4 = x- 2*x.^2+x.^3;

%Convert basis functions to span each interval instead of 0 to 1.
phi_right=subs(f1,x,x/intervalLength);
phi_left=subs(f2,x,x/intervalLength);

psi_right=subs(f3,x,x/intervalLength);
psi_left= subs(f4,x,x/intervalLength);

%Shift the left-handed phi and psi functions 1 interval to the left
phi_left=subs(phi_left,x,x+intervalLength);
psi_left=subs(psi_left,x,x+intervalLength);

%NOTE:At this point phi_left and phi_right are valid on the interval -0.2
%to 0
%phi_right and psi_right are valid on the interval 0 to 0.2

%x_axis = linspace(-1,1,100);
%plot(x_axis,subs(phi_right,x,x_axis),x_axis,subs(phi_left,x,x_axis),x_axis,subs(psi_right,x,x_axis),x_axis,subs(psi_left,x,x_axis))
%axis([-0.2,0.2,-1,1])
%Legend('phi right','phi left', 'psi right', 'psi left')


%Set up empty matrix for holding coefficients
%Put Phi basis functions on odd rows, put Psi basis functions on even rows
matrix1 = zeros(2*(numIntervals+1));

%Create a list of coefficient variables for phi and psi functions
Coeffs = sym(zeros(1, 2*(numIntervals+1)));
for k=1:numIntervals+1
    Coeffs(2*k-1) = sym(sprintf('c%d', k));
    Coeffs(2*k) = sym(sprintf('d%d', k));
end

%Fill center elements of array
%Let first term be summation basis function term, let second term be test
%function

for k = 2:numIntervals
   
    %%Calculate for interval's Phi basis function
    %Phi(k) Phi(k-1)
    matrix1(2*k-1,2*k-3) = int(subs(phi_left,x,x-(k)*intervalLength)*subs(phi_right,x,x-(k-1)*intervalLength),(k-1)*intervalLength, k*intervalLength);
    
    %Phi(k) Psi(k-1)
    matrix1(2*k-1,2*k-2) = int(subs(phi_left,x,x-(k)*intervalLength)*subs(psi_right,x,x-(k-1)*intervalLength),(k-1)*intervalLength, k*intervalLength);
   
    %Phi(k) Phi(k)
    matrix1(2*k-1,2*k-1) = int(subs(phi_left,x,x-(k)*intervalLength)*subs(phi_left,x,x-(k)*intervalLength),(k-1)*intervalLength,k*intervalLength)+int(subs(phi_right,x,x-(k)*intervalLength)*subs(phi_right,x,x-(k)*intervalLength),k*intervalLength, (k+1)*intervalLength);
    
    %Phi(k) Psi(k)
    matrix1(2*k-1,2*k) = int(subs(phi_left,x,x-(k)*intervalLength)*subs(psi_left,x,x-(k)*intervalLength),(k-1)*intervalLength,k*intervalLength)+int(subs(phi_right,x,x-(k)*intervalLength)*subs(psi_right,x,x-(k)*intervalLength),k*intervalLength, (k+1)*intervalLength);
    
    %Phi(k) Phi(k+1)
    matrix1(2*k-1,2*k+1) = int(subs(phi_right,x,x-(k)*intervalLength)*subs(phi_left,x,x-(k+1)*intervalLength),k*intervalLength, (k+1)*intervalLength);
    
    %Phi(k) Psi(k+1)
    matrix1(2*k-1,2*k+2) = int(subs(phi_right,x,x-(k)*intervalLength)*subs(psi_left,x,x-(k+1)*intervalLength),k*intervalLength, (k+1)*intervalLength);   
     
    
    %%Calculate for interval's Psi basis function
    %%Psi(k) Phi(k-1)
    matrix1(2*k,2*k-3) = int(subs(psi_left,x,x-(k)*intervalLength)*subs(phi_right,x,x-(k-1)*intervalLength), (k-1)*intervalLength, k*intervalLength);

    %Psi(k) Psi(k-1)
    matrix1(2*k,2*k-2) = int(subs(psi_left,x,x-(k)*intervalLength)*subs(psi_right,x,x-(k-1)*intervalLength), (k-1)*intervalLength, k*intervalLength);
   
    %Psi(k) Phi(k)
    matrix1(2*k,2*k-1) = int(subs(psi_left,x,x-(k)*intervalLength)*subs(phi_left,x,x-(k)*intervalLength),(k-1)*intervalLength, k*intervalLength)+int(subs(psi_right,x,x-(k)*intervalLength)*subs(phi_right,x,x-(k)*intervalLength),k*intervalLength, (k+1)*intervalLength);
    
    %Psi(k) Psi(k)
    matrix1(2*k,2*k) = int(subs(psi_left,x,x-(k)*intervalLength)*subs(psi_left,x,x-(k)*intervalLength),(k-1)*intervalLength, k*intervalLength)+int(subs(psi_right,x,x-(k)*intervalLength)*subs(psi_right,x,x-(k)*intervalLength),k*intervalLength,(k+1)*intervalLength);
    
    %Psi(k) Phi(k+1)
    matrix1(2*k,2*k+1) = int(subs(psi_right,x,x-(k)*intervalLength)*subs(phi_left,x,x-(k+1)*intervalLength),k*intervalLength, (k+1)*intervalLength);
    
    %Psi(k) Psi(k+1)
    matrix1(2*k,2*k+2) = int(subs(psi_right,x,x-(k)*intervalLength)*subs(psi_left,x,x-(k+1)*intervalLength),k*intervalLength, (k+1)*intervalLength);
    
    
end
%Handle 1st interval 'half case'
%Phis
matrix1(1,1) = int(subs(phi_right,x,x-(0)*intervalLength)*subs(phi_right,x,x-(0)*intervalLength), 0, intervalLength);
matrix1(1,2) = int(subs(phi_right,x,x-(0)*intervalLength)*subs(psi_right,x,x-(0)*intervalLength), 0, intervalLength);
matrix1(1,3) = int(subs(phi_right,x,x-(0)*intervalLength)*subs(phi_left,x,x-(1)*intervalLength),intervalLength, 2*intervalLength);
matrix1(1,4) = int(subs(phi_right,x,x-(0)*intervalLength)*subs(psi_left,x,x-(1)*intervalLength),intervalLength, 2*intervalLength);
    
%Psis
matrix1(2,1) = int(subs(psi_right,x,x-(0)*intervalLength)*subs(phi_right,x,x-(0)*intervalLength), 0, intervalLength);
matrix1(2,2) = int(subs(psi_right,x,x-(0)*intervalLength)*subs(psi_right,x,x-(0)*intervalLength), 0, intervalLength);
matrix1(2,3) = int(subs(psi_right,x,x-(0)*intervalLength)*subs(phi_left,x,x-(1)*intervalLength),intervalLength, 2*intervalLength);
matrix1(2,4) = int(subs(psi_right,x,x-(0)*intervalLength)*subs(psi_left,x,x-(1)*intervalLength),intervalLength, 2*intervalLength);
    
%Handle last interval 'half case'
%Phis
matrix1(2*(numIntervals+1)-1,2*(numIntervals+1)-3)= int(subs(phi_left,x,x-(numIntervals+1)*intervalLength)*subs(phi_right,x,x-(numIntervals)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
matrix1(2*(numIntervals+1)-1,2*(numIntervals+1)-2)= int(subs(phi_left,x,x-(numIntervals+1)*intervalLength)*subs(psi_right,x,x-(numIntervals)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
matrix1(2*(numIntervals+1)-1,2*(numIntervals+1)-1)= int(subs(phi_left,x,x-(numIntervals+1)*intervalLength)*subs(phi_left,x,x-(numIntervals+1)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
matrix1(2*(numIntervals+1)-1,2*(numIntervals+1))= int(subs(phi_left,x,x-(numIntervals+1)*intervalLength)*subs(psi_left,x,x-(numIntervals+1)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
    

%Psis
matrix1(2*(numIntervals+1),2*(numIntervals+1)-3) = int(subs(psi_left,x,x-(numIntervals+1)*intervalLength)*subs(phi_right,x,x-(numIntervals)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
matrix1(2*(numIntervals+1),2*(numIntervals+1)-2) = int(subs(psi_left,x,x-(numIntervals+1)*intervalLength)*subs(psi_right,x,x-(numIntervals)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
matrix1(2*(numIntervals+1),2*(numIntervals+1)-1) = int(subs(psi_left,x,x-(numIntervals+1)*intervalLength)*subs(phi_left,x,x-(numIntervals+1)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
matrix1(2*(numIntervals+1),2*(numIntervals+1)) = int(subs(psi_left,x,x-(numIntervals+1)*intervalLength)*subs(psi_left,x,x-(numIntervals+1)*intervalLength),(numIntervals)*intervalLength, (numIntervals+1)*intervalLength);
     
%Multiply coefficients into LHS matrix
LHS=matrix1

RHS = 12*x^2;
%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(1,2*(numIntervals+1));

for k = 3:2:2*numIntervals
    %Phi test function
    rhs_matrix(k) = int(RHS*subs(phi_left,x,x-(k-1)*intervalLength),k*intervalLength-intervalLength,k*intervalLength)+int(RHS*subs(phi_right,x,x-(k)*intervalLength),k*intervalLength,k*intervalLength+intervalLength);
    
    %Psi test function
    rhs_matrix(k+1) = int(RHS*subs(psi_left,x,x-(k-1)*intervalLength),k*intervalLength-intervalLength,k*intervalLength)+int(RHS*subs(psi_right,x,x-(k)*intervalLength),k*intervalLength,k*intervalLength+intervalLength);
    
    
end


%Handle 1st interval 'half case'
%Phi test function
rhs_matrix(1) = int(RHS*subs(phi_right,x,x-(0)*intervalLength),0,intervalLength);
    
%Psi test function
rhs_matrix(2) = int(RHS*subs(psi_right,x,x-(0)*intervalLength),0,intervalLength);
    

%Handle last interval 'half case'
%Phi test function
rhs_matrix(2*(numIntervals+1)-1) = int(RHS*subs(phi_left,x,x-(numIntervals+1)*intervalLength),(numIntervals)*intervalLength,(numIntervals+1)*intervalLength);
    
%Psi test function
rhs_matrix(2*(numIntervals+1)) = int(RHS*subs(psi_left,x,x-(numIntervals+1)*intervalLength),(numIntervals)*intervalLength,(numIntervals+1)*intervalLength);
    

rhs_matrix

%Solve simultaneous equations for basis function coefficients
Coeffs = mrdivide(rhs_matrix,LHS);


solution = 0;
for index = 3:2:2*numIntervals
    
    solution = solution + Coeffs(index)*subs(phi_left,x,x-(index-1)*intervalLength) + Coeffs(index)*subs(phi_right,x,x-index*intervalLength);
    
    solution = solution + Coeffs(index+1)*subs(psi_left,x,x-(index-1)*intervalLength) + Coeffs(index)*subs(psi_right,x,x-index*intervalLength);
    
end

%Add first and last cases
solution = solution + Coeffs(1)*phi_right + Coeffs(2)*psi_right;

solution = solution + Coeffs(2*numIntervals+1)*subs(phi_left,x,x-(index-1)*intervalLength) + Coeffs(2*numIntervals+2)*subs(psi_left,x,x-(index-1)*intervalLength);



y=linspace(0,1,100);

plot(y,subs(solution,x,y),y,y.^4)

