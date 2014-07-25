clear all
close all
t=cputime;

%Version 1.0 of actual calculations
%Calculates the approximate solution of a known test function
%Uses numerical integration instead of symbolic integration
%The length of the testing regions can be varied
%The number of intervals can be varied
%The starting point of the testing region can be varied

%Collection of constants which may be useful

chi = 4.07; %(*electron affinity in eV*)
gamma = 5.47;  %(*ionization energy in eV*)
phi = chi + gamma/2; %(*work function in eV*)
Eg = gamma - chi;
ni = 2*10^6; %(*intrinsic carrier concentration in [1/cm^3]*)
me = 9.1*10^-31; %(*mass of electron  in kg*)
meffe = .067*me ;%(*mass in density of states for electrons [1]*)
meffh = .47*me; %(*mass in density of states for holes [1]*)
keV = 8.617*10^-5 ;%(*boltzmann  constant in eV/K*)
ksi = 1.38*10^-23; %(*boltzmann constant in J/K*)
T = 300; %(*Cell temperature in K*)
h = 6.626 * 10^-34;
Ec = -chi;
Ev = -gamma;
e0m = 8.854*10^-12;  %(*permittivity of free space in [s^4 A^2/m^3 kg]*)
e0 = e0m*10^-6; %(*permittivity of free space in [s^4 A^2/cm^3 kg]*)
er = 13.1; %(*permittivity ratio for a given material [1]*)
q = 1.6*10^-19; %(*Fundamental charge in C*)
taun = 10^-9; %(*electron lifetime in s*)
taup = 10^-9;%(*hole lifetime in s*)
mun = 8500; %(*electron mobility in cm^2/V s*)
mup = 400; %(*hole mobility in cm^2/V s*)
Dn = 220; %(*electron diffusivity in cm^2/s*)
Dp = 10; %(*hole diffusivity in cm^2/s*)
Sr = 10^5; %(*recombination velocity (in cm/s)*)
Wn = 10^-4; %(*Lenght of n-type region (in cm)*)
Wp = 10^-4; %(*Length of p-type region (in cm)*)


%Input the start and end points of the testing region and determine its
%length
startPoint = 0;
endPoint = Wn;
totalLength = endPoint-startPoint;

%Set number of intervals to divide domain into
numIntervals = 20;
intervalLength = totalLength/numIntervals;

syms x

%(*Solar incidence data*)
alphaOpt = 10000; %(*/cm*) (*GaAs*)
Conc = 1;
Lflux = 2.0*10^17; %(*/cm^2 s *) (*photons absorbed by GaAs at 1 sun*)
R = 0.1; %(*reflection coefficient*)
Gop = alphaOpt*Conc*Lflux*(1 - R)*exp(-alphaOpt*(x));%(*/cm^3*) 


%Set RHS equation to forcing function
RHS = -Gop;


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

%Generate basis functions for each node (index-1), appropriately scaled and
%offset
for index=1:numIntervals+1
phi_right(index)= subs(f1,x,(x-(index-1)*intervalLength-startPoint)/intervalLength);
phi_left(index)=subs(f2,x,(x-(index-2)*intervalLength-startPoint)/intervalLength);

psi_right(index)=subs(f4,x,(x-(index-1)*intervalLength-startPoint)/intervalLength)*intervalLength;
psi_left(index)=subs(f3,x,(x-(index-2)*intervalLength-startPoint)/intervalLength)*intervalLength;
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
%     testPsL(1,index) = subs(psi_left(index),x,(index-2)*intervalLength+startPoint);
%     testPsL(2,index) = subs(psi_left(index),x,(index-1)*intervalLength+startPoint);
%     testPsL(3,index) = subs(diff(psi_left(index),x),x,(index-2)*intervalLength+startPoint);
%     testPsL(4,index) = subs(diff(psi_left(index),x),x,(index-1)*intervalLength+startPoint);
%     
%     
% end

%  x_axis = linspace(-totalLength,totalLength,1000);
%  testIndex=1;
%  plot(x_axis,subs(phi_right(testIndex),x,x_axis),x_axis,subs(phi_left(testIndex),x,x_axis),x_axis,subs(psi_right(testIndex),x,x_axis),x_axis,subs(psi_left(testIndex),x,x_axis))
%  axis([-intervalLength,intervalLength,-1,1])
%  figure
%  plot(x_axis,subs(psi_right(testIndex),x,x_axis),x_axis,subs(psi_left(testIndex),x,x_axis))
%  axis([-intervalLength,intervalLength,-intervalLength,intervalLength])
%  figure



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert initial guess for p(x) to discrete form to be compliant with
%format produced by iterative loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Maximum tolerable error in numeric integration steps, can be this large
%due to magnitude of most outputs (10^(large number))
Error = 100;

%Set up initial guess for p(x).  The curve goes has the vague shape of the
%expected output (peak near 0 and vaguely exponential).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: If pGuess is change, initial BCs for pGuess must also be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1=-2.5*10^11;
p2=1.35*10^12;
p3=-10^8;
gammap1=100000;
gammap2=10000;
gammap3=1;

pGuess = p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(gammap3);

%Convert initial guess into form finite element composition
%This will allow it to be usable in later calculations

%Boundary conditions of test function put into form usable for the
%calculations below

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alpha = subs(diff(pGuess,x),x,0)/subs(pGuess,x,0);

%beta is boundary condition term, n'(L) = beta n(L)
beta= subs(diff(pGuess,x),x,Wn)/subs(pGuess,x,Wn);

matrix = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);
 
%Create a list of coefficient variables for phi and psi functions

%Fill center elements of array
%Let first term be summation basis function term, let second term be test
%function

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix(2*k-4,2*k-2) = quad(matlabFunction(phi_left(k)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Phi(k) Psi(k-1)
    matrix(2*k-3,2*k-2) = quad(matlabFunction(phi_left(k)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Phi(k) Phi(k)
    matrix(2*k-2,2*k-2) = quad(matlabFunction(phi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(phi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k)
    matrix(2*k-1,2*k-2) = quad(matlabFunction(phi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(phi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Phi(k+1)
    matrix(2*k,2*k-2) = quad(matlabFunction(phi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k+1)
    matrix(2*k+1,2*k-2) = quad(matlabFunction(phi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix(2*k-4,2*k-1) = quad(matlabFunction(psi_left(k)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Psi(k) Psi(k-1)
    matrix(2*k-3,2*k-1) = quad(matlabFunction(psi_left(k)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Psi(k) Phi(k)
    matrix(2*k-2,2*k-1) = quad(matlabFunction(psi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(psi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k)
    matrix(2*k-1,2*k-1) = quad(matlabFunction(psi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(psi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi(k) Phi(k+1)
    matrix(2*k,2*k-1) = quad(matlabFunction(psi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k+1)
    matrix(2*k+1,2*k-1) = quad(matlabFunction(psi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
end

%Handle 1st interval 'half case'
%Phis
matrix(1,1) = quad(matlabFunction(phi_right(1)*(phi_right(1)+alpha*psi_right(1))), startPoint, intervalLength+startPoint)+alpha*quad(matlabFunction(psi_right(1)*(phi_right(1)+alpha*psi_right(1))),startPoint, intervalLength+startPoint);
matrix(2,1) = quad(matlabFunction(phi_right(1)*phi_left(2)),startPoint, intervalLength+startPoint)+alpha*quad(matlabFunction(psi_right(1)*phi_left(2)),startPoint, intervalLength+startPoint);
matrix(3,1) = quad(matlabFunction(phi_right(1)*psi_left(2)),startPoint, intervalLength+startPoint)+alpha*quad(matlabFunction(psi_right(1)*psi_left(2)),startPoint, intervalLength+startPoint);

%Handle last interval 'half case'
%Phis
matrix(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(matlabFunction(phi_left(numIntervals+1)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(matlabFunction(psi_left(numIntervals+1)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrix(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(matlabFunction(phi_left(numIntervals+1)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(matlabFunction(psi_left(numIntervals+1)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrix(2*(numIntervals),2*(numIntervals+1)-2)= quad(matlabFunction(phi_left(numIntervals+1)*(phi_left(numIntervals+1)+beta*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(matlabFunction(psi_left(numIntervals+1)*(phi_left(numIntervals+1)+beta*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);

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
matrix(2*k-3,2*k-2) = quad(matlabFunction(phi_left(k)*(phi_right(1)+alpha*psi_right(1))),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix(2*k-2,2*k-2) = quad(matlabFunction(phi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(phi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrix(2*k-1,2*k-2) = quad(matlabFunction(phi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(phi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrix(2*k,2*k-2) = quad(matlabFunction(phi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k+1)
matrix(2*k+1,2*k-2) = quad(matlabFunction(phi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix(2*k-3,2*k-1) = quad(matlabFunction(psi_left(k)*(phi_right(1)+alpha*psi_right(1))), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
  
%Psi(k) Phi(k)
matrix(2*k-2,2*k-1) = quad(matlabFunction(psi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrix(2*k-1,2*k-1) = quad(matlabFunction(psi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(psi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrix(2*k,2*k-1) = quad(matlabFunction(psi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Psi(k) Psi(k+1)
matrix(2*k+1,2*k-1) = quad(matlabFunction(psi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix(2*k-4,2*k-2) = quad(matlabFunction(phi_left(k)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1) 
matrix(2*k-3,2*k-2) = quad(matlabFunction(phi_left(k)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Phi(k)
matrix(2*k-2,2*k-2) = quad(matlabFunction(phi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(phi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrix(2*k-1,2*k-2) = quad(matlabFunction(phi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(phi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrix(2*k,2*k-2) = quad(matlabFunction(phi_right(k)*(phi_left(numIntervals+1)+beta*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix(2*k-4,2*k-1) = quad(matlabFunction(psi_left(k)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

%Psi(k) Psi(k-1)
matrix(2*k-3,2*k-1) = quad(matlabFunction(psi_left(k)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
%Psi(k) Phi(k)
matrix(2*k-2,2*k-1) = quad(matlabFunction(psi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(psi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrix(2*k-1,2*k-1) = quad(matlabFunction(psi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(psi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrix(2*k,2*k-1) = quad(matlabFunction(psi_right(k)*(phi_left(numIntervals+1)+beta*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    

%Multiply coefficients into LHS matrix
LHS=matrix;

%Calculates for RHS of equation (forcing function time test functions)

%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrix(2*k-2) = quad(matlabFunction(pGuess*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(pGuess*phi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi test function
    rhs_matrix(2*k-1) = quad(matlabFunction(pGuess*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction(pGuess*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    
end

%Handle 1st interval 'half case'
%Phi test function
rhs_matrix(1) = quad(matlabFunction(pGuess*(phi_right(1)+alpha*psi_right(1))),startPoint,intervalLength+startPoint);    

%Handle last interval 'half case'
%Phi test function
rhs_matrix(2*(numIntervals+1)-2) = quad(matlabFunction(pGuess*(phi_left(numIntervals+1)+beta*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint,(numIntervals)*intervalLength+startPoint);
   
RHSOut =rhs_matrix;

%Solve simultaneous equations for basis function coefficients
Coeffs = LHS\RHSOut;


%Collect Cj's in an array for display purposes
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);

%Cpjs(1,:)= linspace(startPoint,totalLength+startPoint,numIntervals+1);
%Dpjs(1,:)= linspace(startPoint,totalLength+startPoint,numIntervals+1);
for index=2:length(Cpjs)-1
   Cpjs(1,index) = Coeffs(2*index-2);
   Dpjs(1,index) = Coeffs(2*index-1);
end

Cpjs(1,1) = Coeffs(1);
Dpjs(1,1) = alpha*Coeffs(1);
Cpjs(1,numIntervals+1) = Coeffs(length(Coeffs));
Dpjs(1,numIntervals+1) = beta*Coeffs(length(Coeffs));

testc = Cpjs;
testd = Dpjs;
plot(Cpjs)
figure
%%Cpjs and Dpjs now contain coefficients are each node for p(x) and p'(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End finding initial p(x) finite element coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General BCs for calculated solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alphan = Sr/Dn;

%beta is boundary condition term, n'(L) = beta n(L)
betan = -Sr/Dn;

%alpha is boundary conditions term, p'(0)=alpha*p(0)
alphap = Sr/Dp;
%beta is boundary condition term, p'(L) = beta p(L)
betap = -Sr/Dp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates for RHS of equation (forcing function times test functions)
%This calculation is performed outside the iterative loop, as it will not
%change over the course of changing p(x), n(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrix(2*k-2) = quad(matlabFunction(RHS*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(RHS*phi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi test function
    rhs_matrix(2*k-1) = quad(matlabFunction(RHS*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(RHS*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    
end

%Handle 1st interval 'half case'
%Phi test function
rhs_matrix(1) = quad(matlabFunction(RHS*(phi_right(1)+alphan*psi_right(1))),startPoint,intervalLength+startPoint, Error);    

%Handle last interval 'half case'
%Phi test function
rhs_matrix(2*(numIntervals+1)-2) = quad(matlabFunction(RHS*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint,(numIntervals)*intervalLength+startPoint,Error);
   
RHSOut =rhs_matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End RHS calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with n''(x) term
%Performed outside iterative loop because this term is not influenced by
%p(x) or n(x).  Should not change over calculations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Front end coefficient of n''(x) term
Cn1=Dn - ksi*T*mun/q;

%Set up empty matrix for holding coefficients
%Put Phi basis functions on odd rows, put Psi basis functions on even rows
%Subtract 2 rows and 2 columns to account for Boundary Condition
%limitations
matrix1 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix1(2*k-4,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix1(2*k-3,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k)
    matrix1(2*k-2,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(phi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix1(2*k-1,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(phi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix1(2*k,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix1(2*k+1,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix1(2*k-4,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix1(2*k-3,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix1(2*k-2,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(psi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix1(2*k-1,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(psi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix1(2*k,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix1(2*k+1,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis
matrix1(1,1) = quad(matlabFunction(diff(phi_right(1),x,2)*(phi_right(1)+alphan*psi_right(1))), startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction(diff(psi_right(1),x,2)*(phi_right(1)+alphan*psi_right(1))),startPoint, intervalLength+startPoint,Error);
matrix1(2,1) = quad(matlabFunction(diff(phi_right(1),x,2)*phi_left(2)),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction(diff(psi_right(1),x,2)*phi_left(2)),startPoint, intervalLength+startPoint,Error);
matrix1(3,1) = quad(matlabFunction(diff(phi_right(1),x,2)*psi_left(2)),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction(diff(psi_right(1),x,2)*psi_left(2)),startPoint, intervalLength+startPoint,Error);

quad(matlabFunction(diff(phi_right(1),x,2)*phi_right(1)), startPoint, intervalLength+startPoint,Error)
alphan*quad(matlabFunction(diff(psi_right(1),x,2)*phi_right(1)),startPoint, intervalLength+startPoint,Error)

betan
%Handle last interval 'half case'
%Phis
matrix1(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(matlabFunction(diff(phi_left(numIntervals+1),x,2)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction(diff(psi_left(numIntervals+1),x,2)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix1(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(matlabFunction(diff(phi_left(numIntervals+1),x,2)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction(diff(psi_left(numIntervals+1),x,2)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix1(2*(numIntervals),2*(numIntervals+1)-2)= quad(matlabFunction(diff(phi_left(numIntervals+1),x,2)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction(diff(psi_left(numIntervals+1),x,2)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);

quad(matlabFunction(diff(phi_left(numIntervals+1),x,2)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)
betan*quad(matlabFunction(diff(psi_left(numIntervals+1),x,2)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)


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
matrix1(2*k-3,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*(phi_right(1)+alphan*psi_right(1))),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(phi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(phi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix1(2*k+1,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-3,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*(phi_right(1)+alphan*psi_right(1))), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(psi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(psi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix1(2*k+1,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix1(2*k-4,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix1(2*k-3,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(phi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(matlabFunction(diff(phi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(phi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) = quad(matlabFunction(diff(phi_right(k),x,2)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-4,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix1(2*k-3,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(psi_right(k),x,2)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(matlabFunction(diff(psi_left(k),x,2)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(psi_right(k),x,2)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) = quad(matlabFunction(diff(psi_right(k),x,2)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

Matrix1Final = Cn1*matrix1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of n''(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with n'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Front end coefficient for this matrix.  DOES NOT INCLUDE CONTRIBUTION FROM
%p(x) DISTRIBUTION.
Cn2 = ksi*T*mun/q;

matrix2 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%%pterm = (p'(x)/p(x))
%%term associated with phi_left(k) interval of n(x)
%%diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)
%%term associated with phi_right(k) interval of n(x)
%%diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)
for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Cpk
    
    %Phi(k) Phi(k-1)    
    matrix2(2*k-4,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix2(2*k-3,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k)
    matrix2(2*k-2,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix2(2*k-1,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix2(2*k,2*k-2) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix2(2*k+1,2*k-2) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix2(2*k-4,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix2(2*k-3,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix2(2*k-2,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix2(2*k-1,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix2(2*k,2*k-1) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix2(2*k+1,2*k-1) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis
matrix2(1,1) = quad(matlabFunction(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(1),x,1)*(phi_right(1)+alphan*psi_right(1))),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(1),x,1)*(phi_right(1)+alphan*psi_right(1))),startPoint, intervalLength+startPoint,Error);
matrix2(2,1) = quad(matlabFunction(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(1),x,1)*phi_left(2)),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(1),x,1)*phi_left(2)),startPoint, intervalLength+startPoint,Error);
matrix2(3,1) = quad(matlabFunction(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(1),x,1)*psi_left(2)),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(1),x,1)*psi_left(2)),startPoint, intervalLength+startPoint,Error);


%Handle last interval 'half case'
%Phis
matrix2(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(matlabFunction(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(numIntervals+1),x,1)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(numIntervals+1),x,1)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix2(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(matlabFunction(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(numIntervals+1),x,1)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(numIntervals+1),x,1)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix2(2*(numIntervals),2*(numIntervals+1)-2)= quad(matlabFunction(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(numIntervals+1),x,1)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(numIntervals+1),x,1)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);

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
matrix2(2*k-3,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*(phi_right(1)+alphan*psi_right(1))),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix2(2*k-2,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix2(2*k-1,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix2(2*k,2*k-2) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix2(2*k+1,2*k-2) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix2(2*k-3,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*(phi_right(1)+alphan*psi_right(1))), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix2(2*k-2,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix2(2*k-1,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix2(2*k,2*k-1) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix2(2*k+1,2*k-1) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   

%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix2(2*k-4,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix2(2*k-3,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix2(2*k-2,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix2(2*k-1,2*k-2) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_left(k),x,1)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix2(2*k,2*k-2) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(phi_right(k),x,1)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix2(2*k-4,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix2(2*k-3,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix2(2*k-2,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix2(2*k-1,2*k-1) = quad(matlabFunction(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_left(k),x,1)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix2(2*k,2*k-1) = quad(matlabFunction(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)*diff(psi_right(k),x,1)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);

Matrix2Final = Cn2*matrix2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of n'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with n(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix3 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%Coefficient Cn3 is built into filling the matrix since it has a difference
%that can't be removed built-in.

%%pterm = ksi*T*mun/q*(p''(x)/p(x)-(p'(x)/p(x))^2)-1/taun
%%term associated with phi_left(k) interval of n(x)
%%(ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)
%%term associated with phi_right(k) interval of n(x)
%%(ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix3(2*k-4,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix3(2*k-3,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k)
    matrix3(2*k-2,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix3(2*k-1,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Phi(k+1)
    matrix3(2*k,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix3(2*k+1,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix3(2*k-4,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix3(2*k-3,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix3(2*k-2,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix3(2*k-1,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix3(2*k,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix3(2*k+1,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis
matrix3(1,1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(1)*(phi_right(1)+alphan*psi_right(1))),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction((ksi*T*mun/q*(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(1)*(phi_right(1)+alphan*psi_right(1))),startPoint, intervalLength+startPoint,Error);
matrix3(2,1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(1)*phi_left(2)),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction((ksi*T*mun/q*(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(1)*phi_left(2)),startPoint, intervalLength+startPoint,Error);
matrix3(3,1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(1)*psi_left(2)),startPoint, intervalLength+startPoint,Error)+alphan*quad(matlabFunction((ksi*T*mun/q*(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(2,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(1)*psi_left(2)),startPoint, intervalLength+startPoint,Error);

%Handle last interval 'half case'
%Phis
matrix3(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(matlabFunction((ksi*T*mun/q*(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(numIntervals+1)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction((ksi*T*mun/q*(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(numIntervals+1)*phi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix3(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(matlabFunction((ksi*T*mun/q*(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(numIntervals+1)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction((ksi*T*mun/q*(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(numIntervals+1)*psi_right(numIntervals)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix3(2*(numIntervals),2*(numIntervals+1)-2)= quad(matlabFunction((ksi*T*mun/q*(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(numIntervals+1)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(matlabFunction((ksi*T*mun/q*(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(numIntervals+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(numIntervals+1)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);

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
matrix3(2*k-3,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*(phi_right(1)+alphan*psi_right(1))),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix3(2*k-2,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix3(2*k-1,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix3(2*k,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix3(2*k+1,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix3(2*k-3,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*(phi_right(1)+alphan*psi_right(1))), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix3(2*k-2,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix3(2*k-1,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix3(2*k,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*phi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix3(2*k+1,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*psi_left(k+1)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix3(2*k-4,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*phi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
matrix3(2*k-3,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*psi_right(k-1)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix3(2*k-2,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix3(2*k-1,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix3(2*k,2*k-2) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*phi_right(k)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix3(2*k-4,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*phi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix3(2*k-3,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*psi_right(k-1)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix3(2*k-2,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*phi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*phi_right(k)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix3(2*k-1,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_left(k)*psi_left(k)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*psi_right(k)),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix3(2*k,2*k-1) = quad(matlabFunction((ksi*T*mun/q*(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,2)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left) -(diff(px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left),x,1)/px(k+1,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left))^2) - 1/taun)*psi_right(k)*(phi_left(numIntervals+1)+betan*psi_left(numIntervals+1))),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);


Matrix3Final = matrix3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of n(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Combine expressions from each portion for total matrix of LHS
LHS = Matrix1Final +Matrix2Final +Matrix3Final;



%Solve simultaneous equations for basis function coefficients
Coeffs = LHS\RHSOut;

Coeffs1 = Matrix1Final\RHSOut;
Coeffs2 = Matrix2Final\RHSOut;
Coeffs3 = Matrix3Final\RHSOut;

%Collect Cj's in an array for display purposes
guesstimate = zeros(2,numIntervals+1);
guesstimate(1,:)= linspace(startPoint,totalLength+startPoint,numIntervals+1);
for index=2:length(guesstimate)
   guesstimate(2,index) = Coeffs(2*index-2);
end
guesstimate(2,1) = Coeffs(1);

guesstimate1 = zeros(2,numIntervals+1);
guesstimate1(1,:)= linspace(startPoint,totalLength+startPoint,numIntervals+1);
for index=2:length(guesstimate1)
   guesstimate1(2,index) = Coeffs1(2*index-2);
end
guesstimate1(2,1) = Coeffs1(1);


guesstimate2 = zeros(2,numIntervals+1);
guesstimate2(1,:)= linspace(startPoint,totalLength+startPoint,numIntervals+1);
for index=2:length(guesstimate2)
   guesstimate2(2,index) = Coeffs2(2*index-2);
end
guesstimate2(2,1) = Coeffs2(1);

guesstimate3 = zeros(2,numIntervals+1);
guesstimate3(1,:)= linspace(startPoint,totalLength+startPoint,numIntervals+1);
for index=2:length(guesstimate3)
   guesstimate3(2,index) = Coeffs3(2*index-2);
end
guesstimate3(2,1) = Coeffs3(1);

%Plotting
y=linspace(startPoint,totalLength+startPoint,100);

plot(y,subs(pGuess,x,y))

figure
plot(guesstimate(1,:),guesstimate(2,:));title('all Ms');

figure
plot(guesstimate1(1,:),guesstimate1(2,:));title('M1');

figure
plot(guesstimate2(1,:),guesstimate2(2,:));title('M2');

figure
plot(guesstimate3(1,:),guesstimate3(2,:));title('M3');


%Displays runtime of m-file
time = cputime-t