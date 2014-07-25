clear all
close all
t=cputime;
%profile on

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
numIntervals = 100;
intervalLength = totalLength/numIntervals;
maxIterations=1;



%Random noise bounds
randLow=1;
randHigh=1;


filepath=strcat('stability test v4/rand',num2str(randLow),'/',num2str(numIntervals),'terms');
mkdir(filepath)

syms x

%(*Solar incidence data*)
alphaOpt = 10000; %(*/cm*) (*GaAs*)
Conc = 1;
Lflux = 2.0*10^17; %(*/cm^2 s *) (*photons absorbed by GaAs at 1 sun*)
R = 0.1; %(*reflection coefficient*)
Gop = alphaOpt*Conc*Lflux*(1 - R)*exp(-alphaOpt*(x));%(*/cm^3*) 


indexp = linspace(startPoint,totalLength+startPoint,numIntervals+1);
indexn = linspace(startPoint,totalLength+startPoint,numIntervals+1);

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
phi_right=cell(1,numIntervals+1);
phi_left=cell(1,numIntervals+1);
psi_right=cell(1,numIntervals+1);
psi_left=cell(1,numIntervals+1);

%Generate basis functions for each node (index-1), appropriately scaled and
%offset
%Use example call of phi_right{1}(x) to get function in 1st basis function
%in terms of variable x.

for index=1:numIntervals+1
phi_right{index}=matlabFunction(subs(f1,x,(x-(index-1)*intervalLength-startPoint)/intervalLength));
phi_left{index}=matlabFunction(subs(f2,x,(x-(index-2)*intervalLength-startPoint)/intervalLength));

psi_right{index}=matlabFunction(subs(f4,x,(x-(index-1)*intervalLength-startPoint)/intervalLength)*intervalLength);
psi_left{index}=matlabFunction(subs(f3,x,(x-(index-2)*intervalLength-startPoint)/intervalLength)*intervalLength);
end

% x_axis=0:intervalLength/100:intervalLength;
% subplot(221),plot(x_axis,phi_right{1}(x_axis))
% subplot(222),plot(x_axis,phi_left{2}(x_axis))
% subplot(223),plot(x_axis,psi_right{1}(x_axis))
% subplot(224),plot(x_axis,psi_left{2}(x_axis))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Symbolic representation of basis functions
%Required as derivatives of function handles are not possible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create empty arrays to hold basis functions for each node
phi_right_syms=sym(zeros(1,numIntervals+1));
phi_left_syms=sym(zeros(1,numIntervals+1));
psi_right_syms=sym(zeros(1,numIntervals+1));
psi_left_syms=sym(zeros(1,numIntervals+1));

%Generate basis functions for each node (index-1), appropriately scaled and
%offset
for index=1:numIntervals+1
phi_right_syms(index)= subs(f1,x,(x-(index-1)*intervalLength-startPoint)/intervalLength);
phi_left_syms(index)=subs(f2,x,(x-(index-2)*intervalLength-startPoint)/intervalLength);

psi_right_syms(index)=subs(f4,x,(x-(index-1)*intervalLength-startPoint)/intervalLength)*intervalLength;
psi_left_syms(index)=subs(f3,x,(x-(index-2)*intervalLength-startPoint)/intervalLength)*intervalLength;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End symbolic basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute function handles of 1st derivatives of basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_phi_right=cell(1,numIntervals+1);
d_phi_left=cell(1,numIntervals+1);
d_psi_right=cell(1,numIntervals+1);
d_psi_left=cell(1,numIntervals+1);

for index=1:numIntervals+1
d_phi_right{index}= matlabFunction(diff(subs(f1,x,(x-(index-1)*intervalLength-startPoint)/intervalLength),x,1));
d_phi_left{index}=matlabFunction(diff(subs(f2,x,(x-(index-2)*intervalLength-startPoint)/intervalLength),x,1));

d_psi_right{index}=matlabFunction(diff(subs(f4,x,(x-(index-1)*intervalLength-startPoint)/intervalLength)*intervalLength,x,1));
d_psi_left{index}=matlabFunction(diff(subs(f3,x,(x-(index-2)*intervalLength-startPoint)/intervalLength)*intervalLength,x,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute function handles of 2nd derivatives of basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2_phi_right=cell(1,numIntervals+1);
d2_phi_left=cell(1,numIntervals+1);
d2_psi_right=cell(1,numIntervals+1);
d2_psi_left=cell(1,numIntervals+1);

for index=1:numIntervals+1
d2_phi_right{index}= matlabFunction(diff(subs(f1,x,(x-(index-1)*intervalLength-startPoint)/intervalLength),x,2));
d2_phi_left{index}=matlabFunction(diff(subs(f2,x,(x-(index-2)*intervalLength-startPoint)/intervalLength),x,2));

d2_psi_right{index}=matlabFunction(diff(subs(f4,x,(x-(index-1)*intervalLength-startPoint)/intervalLength)*intervalLength,x,2));
d2_psi_left{index}=matlabFunction(diff(subs(f3,x,(x-(index-2)*intervalLength-startPoint)/intervalLength)*intervalLength,x,2));
end


%NOTE:At each node phi_left and psi_left are valid on the interval previous
%to the node, phi_right and psi_right are valid on the interval following
%the node.

% %Test to verify proper performance of basis functions at interval endpoints
% testPhR = zeros(1,numIntervals+1);
% testPhL = zeros(1,numIntervals+1);
% testPsR = zeros(1,numIntervals+1);
% testPsL = zeros(1,numIntervals+1);
% 
% for index = 1:numIntervals+1
%     testPhR(1,index) = subs(phi_right{index},x,(index-1)*intervalLength+startPoint);
%     testPhR(2,index) = subs(phi_right{index},x,(index)*intervalLength+startPoint);
%     testPhR(3,index) = subs(diff(phi_right{index},x),x,(index-1)*intervalLength+startPoint);
%     testPhR(4,index) = subs(diff(phi_right{index},x),x,(index)*intervalLength+startPoint);
%     
%     testPsR(1,index) = subs(psi_right{index},x,(index-1)*intervalLength+startPoint);
%     testPsR(2,index) = subs(psi_right{index},x,(index)*intervalLength+startPoint);
%     testPsR(3,index) = subs(diff(psi_right{index},x),x,(index-1)*intervalLength+startPoint);
%     testPsR(4,index) = subs(diff(psi_right{index},x),x,(index)*intervalLength+startPoint);
%         
%     testPhL(1,index) = subs(phi_left{index},x,(index-2)*intervalLength+startPoint);
%     testPhL(2,index) = subs(phi_left{index},x,(index-1)*intervalLength+startPoint);
%     testPhL(3,index) = subs(diff(phi_left{index},x),x,(index-2)*intervalLength+startPoint);
%     testPhL(4,index) = subs(diff(phi_left{index},x),x,(index-1)*intervalLength+startPoint);
%         
%     testPsL(1,index) = subs(psi_left{index},x,(index-2)*intervalLength+startPoint);
%     testPsL(2,index) = subs(psi_left{index},x,(index-1)*intervalLength+startPoint);
%     testPsL(3,index) = subs(diff(psi_left{index},x),x,(index-2)*intervalLength+startPoint);
%     testPsL(4,index) = subs(diff(psi_left{index},x),x,(index-1)*intervalLength+startPoint);
%     
%     
% end
%
%  x_axis = linspace(-totalLength,totalLength,10000);
%  testIndex=1;
%  plot(x_axis,subs(phi_right{testIndex},x,x_axis),x_axis,subs(phi_left{testIndex},x,x_axis),x_axis,subs(psi_right{testIndex},x,x_axis),x_axis,subs(psi_left{testIndex},x,x_axis))
%  axis([-intervalLength,intervalLength,-1,1])
%  figure
%  plot(x_axis,subs(psi_right{testIndex},x,x_axis),x_axis,subs(psi_left{testIndex},x,x_axis))
%  axis([-intervalLength,intervalLength,-intervalLength,intervalLength])
%  figure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up known n(x) and p(x) distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Maximum tolerable error in numeric integration steps, can be this large
%due to magnitude of most outputs (10^(large number))
Error = 100;

%Set up initial guess for p(x).  The curve goes has the vague shape of the
%expected output (peak near 0 and vaguely exponential).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: If pGuess is change, initial BCs for pGuess must also be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1=-2.4469*10^11;
p2=1.34584*10^12;
p3=-10^9;
gammap1=100000;
gammap2=10000;
gammap3=1;

%Symbolic form of pGuess
pGuess_syms = p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(-gammap3*x);
pprime_syms=diff(pGuess_syms,x);

pGuess = matlabFunction(pGuess_syms);
pprime = matlabFunction(pprime_syms);
pGuee_orig = pGuess;

n1=-1.79577*10^10;
n2=1.29167*10^11;
n3=10^12;
gamman1=100000;
gamman2=10000;
gamman3=1;

%Symbolic form of nGuess
nGuess_syms = n1*exp(-gamman1*x)+n2*exp(-gamman2*x)+n3*exp(-gamman3*x);
nGuess = matlabFunction(nGuess_syms);
nprime_syms = diff(nGuess_syms,x);
nprime = matlabFunction(nprime_syms);


cn1 = Dn - ksi*T*mun/q;
cn2 = ksi*T*mun/q *diff(pGuess_syms, x)/pGuess_syms;
cn3 = ksi*T*mun/q*(diff(pGuess_syms, x,2)/pGuess_syms - (diff(pGuess_syms, x)/pGuess_syms)^2) - 1/taun;

cp1 = ksi*T*mup/q - Dp;
cp2 = -ksi*T*mup/q*diff(nGuess_syms, x)/nGuess_syms;
cp3 = ksi*T*mup/q*((diff(nGuess_syms, x)/nGuess_syms)^2 - diff(nGuess_syms, x,2)/nGuess_syms) - 1/taup;

Gopn_syms = cn1*diff(nGuess_syms,x,2)+cn2*diff(nGuess_syms,x)+cn3*nGuess_syms;
Gopp_syms = cp1*diff(pGuess_syms,x,2)+cp2*diff(pGuess_syms,x)+cp3*pGuess_syms;

Gopn = matlabFunction(Gopn_syms);
Gopp = matlabFunction(Gopp_syms);

%Plotting
y=linspace(startPoint,totalLength+startPoint,100);
fig1 = figure;
plot(y,subs(pGuess_syms,x,y)), title('pGuess')
saveas(fig1,strcat(filepath,'/pGuess.jpg'))

fig2 = figure;
plot(y,subs(nGuess_syms,x,y)), title('nGuess')
saveas(fig2,strcat(filepath,'/nGuess.jpg'))

%Generate random noise instand of pGuess
%Set upper and lower bounds of random number coefficients for each Cj, Dj
p1=-2.4469*10^11*(randLow+(randHigh-randLow)*rand(1));
p2=1.34584*10^12*(randLow+(randHigh-randLow)*rand(1));
p3=-10^9*(randLow+(randHigh-randLow)*rand(1));
gammap1=100000*(randLow+(randHigh-randLow)*rand(1));
gammap2=10000*(randLow+(randHigh-randLow)*rand(1));
gammap3=1*(randLow+(randHigh-randLow)*rand(1));

%Symbolic form of pGuess with noise
%pGuess = p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(-gammap3*x);
offset=0;%Wn/100;
pGuess = matlabFunction(subs(pGuess_syms,x,x-offset));

fig3 = figure;
plot(y,subs(pGuess_syms,x,y),y,subs(pGuess,x,y)), title('pGuess (noise)')
saveas(fig3,strcat(filepath,'/pGuess_noisy.jpg'))

%Convert initial guess into form finite element composition
%This will allow it to be usable in later calculations

%Boundary conditions of test function put into form usable for the
%calculations below

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alpha = subs(diff(pGuess_syms,x),x,0-offset)/subs(pGuess_syms,x,0-offset);

%beta is boundary condition term, n'(L) = beta n(L)
beta= subs(diff(pGuess_syms,x),x,Wn-offset)/subs(pGuess_syms,x,Wn-offset);

matrix = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);
 
%Create a list of coefficient variables for phi and psi functions

%Fill center elements of array
%Let first term be summation basis function term, let second term be test
%function

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix(2*k-4,2*k-2) = quad(@(x) phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Phi(k) Psi(k-1)
    matrix(2*k-3,2*k-2) = quad(@(x) phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Phi(k) Phi(k)
    matrix(2*k-2,2*k-2) = quad(@(x) phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k)
    matrix(2*k-1,2*k-2) = quad(@(x) phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Phi(k+1) 
    matrix(2*k,2*k-2) =   quad(@(x) phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k+1)
    matrix(2*k+1,2*k-2) = quad(@(x) phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix(2*k-4,2*k-1) = quad(@(x) psi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Psi(k) Psi(k-1)
    matrix(2*k-3,2*k-1) = quad(@(x) psi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Psi(k) Phi(k)
    matrix(2*k-2,2*k-1) = quad(@(x) psi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k)
    matrix(2*k-1,2*k-1) = quad(@(x) psi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi(k) Phi(k+1)
    matrix(2*k,2*k-1) =   quad(@(x) psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k+1)
    matrix(2*k+1,2*k-1) = quad(@(x) psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
end

%Handle 1st interval 'half case'
%Phis
matrix(1,1) = quad(@(x) phi_right{1}(x).*(phi_right{1}(x)+alpha.*psi_right{1}(x)), startPoint, intervalLength+startPoint)+alpha*quad(@(x) psi_right{1}(x).*(phi_right{1}(x)+alpha.*psi_right{1}(x)),startPoint, intervalLength+startPoint);
matrix(2,1) = quad(@(x) phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint)                           +alpha*quad(@(x) psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint);
matrix(3,1) = quad(@(x) phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint)                           +alpha*quad(@(x) psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint);

%Handle last interval 'half case'
%Phis 
matrix(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(@(x) psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrix(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(@(x) psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrix(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+beta.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+beta*quad(@(x) psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+beta.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);

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
matrix(2*k-3,2*k-2) = quad(@(x) phi_left{k}(x).*(phi_right{1}(x)+alpha.*psi_right{1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix(2*k-2,2*k-2) = quad(@(x) phi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrix(2*k-1,2*k-2) = quad(@(x) phi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrix(2*k,2*k-2) =   quad(@(x) phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k+1)
matrix(2*k+1,2*k-2) = quad(@(x) phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix(2*k-3,2*k-1) = quad(@(x) psi_left{k}(x).*(phi_right{1}(x)+alpha.*psi_right{1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
  
%Psi(k) Phi(k)
matrix(2*k-2,2*k-1) = quad(@(x) psi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrix(2*k-1,2*k-1) = quad(@(x) psi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrix(2*k,2*k-1) =   quad(@(x) psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Psi(k) Psi(k+1)
matrix(2*k+1,2*k-1) = quad(@(x) psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix(2*k-4,2*k-2) = quad(@(x) phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1) 
matrix(2*k-3,2*k-2) = quad(@(x) phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Phi(k)
matrix(2*k-2,2*k-2) = quad(@(x) phi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrix(2*k-1,2*k-2) = quad(@(x) phi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrix(2*k,2*k-2) =   quad(@(x) phi_right{k}(x).*(phi_left{k+1}(x)+beta.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix(2*k-4,2*k-1) = quad(@(x) psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

%Psi(k) Psi(k-1)
matrix(2*k-3,2*k-1) = quad(@(x) psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
%Psi(k) Phi(k)
matrix(2*k-2,2*k-1) = quad(@(x) psi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrix(2*k-1,2*k-1) = quad(@(x) psi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrix(2*k,2*k-1) =   quad(@(x) psi_right{k}(x).*(phi_left{k+1}(x)+beta.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    

%Multiply coefficients into LHS matrix
LHStest=matrix;

%Calculates for RHS of equation (forcing function time test functions)

%Generate matrix of test functions against RHS of equation
rhs_matrixp = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrixp(2*k-2) = quad(@(x) pGuess(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) pGuess(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi test function
    rhs_matrixp(2*k-1) = quad(@(x) pGuess(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) pGuess(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    
end

%Handle 1st interval 'half case'
%Phi test function 
rhs_matrixp(1) = quad(@(x) pGuess(x).*(phi_right{1}(x)+alpha.*psi_right{1}(x)),startPoint,intervalLength+startPoint);    

%Handle last interval 'half case'
%Phi test function 
rhs_matrixp(2*(numIntervals+1)-2) = quad(@(x) pGuess(x).*(phi_left{numIntervals+1}(x)+beta.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint,(numIntervals)*intervalLength+startPoint);
   
RHSOutp =rhs_matrixp;

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alphatestn = subs(diff(nGuess,x),x,0)/subs(nGuess,x,0);

%beta is boundary condition term, n'(L) = beta n(L)
betatestn= subs(diff(nGuess,x),x,Wn)/subs(nGuess,x,Wn);

matrixn = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);
 
%Create a list of coefficient variables for phi and psi functions

%Fill center elements of array
%Let first term be summation basis function term, let second term be test
%function

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrixn(2*k-4,2*k-2) = quad(@(x) phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Phi(k) Psi(k-1)
    matrixn(2*k-3,2*k-2) = quad(@(x) phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Phi(k) Phi(k)
    matrixn(2*k-2,2*k-2) = quad(@(x) phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k)
    matrixn(2*k-1,2*k-2) = quad(@(x) phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Phi(k+1) 
    matrixn(2*k,2*k-2) =   quad(@(x) phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Phi(k) Psi(k+1)
    matrixn(2*k+1,2*k-2) = quad(@(x) phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrixn(2*k-4,2*k-1) = quad(@(x) psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

    %Psi(k) Psi(k-1)
    matrixn(2*k-3,2*k-1) = quad(@(x) psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
    %Psi(k) Phi(k)
    matrixn(2*k-2,2*k-1) = quad(@(x) psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k)
    matrixn(2*k-1,2*k-1) = quad(@(x) psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi(k) Phi(k+1)
    matrixn(2*k,2*k-1) =   quad(@(x) psi_right{k}(x).*phi_left{k+1}(x), (k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    %Psi(k) Psi(k+1)
    matrixn(2*k+1,2*k-1) = quad(@(x) psi_right{k}(x).*psi_left{k+1}(x), (k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
end

%Handle 1st interval 'half case'
%Phis
matrixn(1,1) = quad(@(x) phi_right{1}(x).*(phi_right{1}(x)+alphatestn.*psi_right{1}(x)), startPoint, intervalLength+startPoint)+alphatestn*quad(@(x) psi_right{1}(x).*(phi_right{1}(x)+alphatestn.*psi_right{1}(x)),startPoint, intervalLength+startPoint);
matrixn(2,1) = quad(@(x) phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint)                                +alphatestn*quad(@(x) psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint);
matrixn(3,1) = quad(@(x) phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint)                                +alphatestn*quad(@(x) psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint);

%Handle last interval 'half case'
%Phis 
matrixn(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)                                          +betatestn*quad(@(x) psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrixn(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)                                          +betatestn*quad(@(x) psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);
matrixn(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betatestn.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint)+betatestn*quad(@(x) psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betatestn.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint);

%Handle interval 2 and n-1 case
%These intervals generate terms in the rows that would normally be in the
%for-loop matrixn, but are omitted due to boundary conditions.  Due to the
%unusual number fo rows included in the final matrixn, they must be input
%separately.

%Interval 2
%Phis
%Code is copied for code in for-loop above
k=2;
%Phi(k) Phi(k-1)
matrixn(2*k-3,2*k-2) = quad(@(x) phi_left{k}(x).*(phi_right{1}(x)+alphatestn.*psi_right{1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrixn(2*k-2,2*k-2) = quad(@(x) phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrixn(2*k-1,2*k-2) = quad(@(x) phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrixn(2*k,2*k-2) =   quad(@(x) phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k+1)
matrixn(2*k+1,2*k-2) = quad(@(x) phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrixn(2*k-3,2*k-1) = quad(@(x) psi_left{k}(x).*(phi_right{1}(x)+alphatestn.*psi_right{1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
  
%Psi(k) Phi(k)
matrixn(2*k-2,2*k-1) = quad(@(x) psi_left{k}(x).* phi_left{k}(x),  (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrixn(2*k-1,2*k-1) = quad(@(x) psi_left{k}(x).* psi_left{k}(x),  (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrixn(2*k,2*k-1) =   quad(@(x) psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Psi(k) Psi(k+1)
matrixn(2*k+1,2*k-1) = quad(@(x) psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrixn(2*k-4,2*k-2) = quad(@(x) phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Psi(k-1) 
matrixn(2*k-3,2*k-2) = quad(@(x) phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
 
%Phi(k) Phi(k)
matrixn(2*k-2,2*k-2) = quad(@(x) phi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Psi(k)
matrixn(2*k-1,2*k-2) = quad(@(x) phi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%Phi(k) Phi(k+1)
matrixn(2*k,2*k-2) =   quad(@(x) phi_right{k}(x).*(phi_left{k+1}(x)+betatestn.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    
%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrixn(2*k-4,2*k-1) = quad(@(x) psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);

%Psi(k) Psi(k-1)
matrixn(2*k-3,2*k-1) = quad(@(x) psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint);
   
%Psi(k) Phi(k)
matrixn(2*k-2,2*k-1) = quad(@(x) psi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
   
%Psi(k) Psi(k)
matrixn(2*k-1,2*k-1) = quad(@(x) psi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint)+quad(@(x) psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
%Psi(k) Phi(k+1)
matrixn(2*k,2*k-1) =   quad(@(x) psi_right{k}(x).*(phi_left{k+1}(x)+betatestn.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint);
    

%Multiply coefficients into LHS matrix
LHStestn=matrixn;

%Generate matrix of test functions against RHS of equation
rhs_matrixn = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrixn(2*k-2) = quad(@(x) nGuess(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) nGuess(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    %Psi test function
    rhs_matrixn(2*k-1) = quad(@(x) nGuess(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint)+quad(@(x) nGuess(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint);
    
    
end

%Handle 1st interval 'half case'
%Phi test function 
rhs_matrixn(1) = quad(@(x) nGuess(x).*(phi_right{1}(x)+alphatestn.*psi_right{1}(x)),startPoint,intervalLength+startPoint);    

%Handle last interval 'half case'
%Phi test function 
rhs_matrixn(2*(numIntervals+1)-2) = quad(@(x) nGuess(x).*(phi_left{numIntervals+1}(x)+betatestn.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint,(numIntervals)*intervalLength+startPoint);
  
RHSOutn =rhs_matrixn;

RHStest=rhs_matrixn;

%Solve simultaneous equations for basis function coefficients
Coeffsp = LHStest\RHSOutp;
Coeffsn = LHStestn\RHSOutn;

%Collect Cj's in an array for display purposes
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);



for index=2:length(Cpjs)-1
   Cpjs(1,index) = Coeffsp(2*index-2);
   Dpjs(1,index) = Coeffsp(2*index-1);
end

Cpjs(1,1) = Coeffsp(1);
Dpjs(1,1) = alpha*Coeffsp(1);


Cpjs(1,numIntervals+1) = Coeffsp(length(Coeffsp));
Dpjs(1,numIntervals+1) = beta*Coeffsp(length(Coeffsp));


%Collect Cj's in an array for display purposes
Cnjs = zeros(1,numIntervals+1);
Dnjs = zeros(1,numIntervals+1);


for index=2:length(Cpjs)-1
   Cnjs(1,index) = Coeffsn(2*index-2);
   Dnjs(1,index) = Coeffsn(2*index-1);
end

Cnjs(1,1) = Coeffsn(1);
Dnjs(1,1) = alphatestn*Coeffsn(1);
Cnjs(1,numIntervals+1) = Coeffsn(length(Coeffsn));
Dnjs(1,numIntervals+1) = betatestn*Coeffsn(length(Coeffsn));

InitialCpjs = Cpjs;
InitialDpjs = Dpjs;

InitialCnjs = Cnjs;
InitialDnjs = Dnjs;

fig3 = figure;
plot(indexp,Cpjs,indexp,subs(pGuess,x,indexp)),title('pGuess (recomposition)');
saveas(fig3,strcat(filepath,'/pGuess_reomposition.jpg'))

fig4 = figure;
plot(indexn,Cnjs,indexn,subs(nGuess,x,indexn)),title('nGuess (recomposition)');
saveas(fig4,strcat(filepath,'/nGuess_reomposition.jpg'))

fig9 = figure;
plot(indexp,Dpjs,indexp,subs(diff(pGuess,x,1),x,indexp)),title('pprime (recomposition)');
saveas(fig9,strcat(filepath,'/pprime_reomposition.jpg'))


fig10 = figure;
plot(indexn,Dnjs,indexn,subs(diff(nGuess,x,1),x,indexn)),title('nprime (recomposition)');
saveas(fig10,strcat(filepath,'/nprime_recomposition.jpg'))

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

pGuess = matlabFunction(pGuess_syms);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates for RHS of equation (forcing function times test functions)
%This calculation is performed outside the iterative loop, as it will not
%change over the course of changing p(x), n(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RHS = Gopn;

%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrix(2*k-2) = quad(@(x) RHS(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) RHS(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi test function
    rhs_matrix(2*k-1) = quad(@(x) RHS(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) RHS(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    
end


%Handle 1st interval 'half case'
%Phi test function
rhs_matrix(1) = quad(@(x) RHS(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),startPoint,intervalLength+startPoint, Error);    

%Handle last interval 'half case'
%Phi test function
rhs_matrix(2*(numIntervals+1)-2) = quad(@(x) RHS(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint,(numIntervals)*intervalLength+startPoint,Error);
   
RHSOutn =rhs_matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End RHS calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with n''(x) term
% %Performed outside iterative loop because this term is not influenced by
% %p(x) or n(x).  Should not change over calculations.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Front end coefficient of n''(x) term
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
    matrix1(2*k-4,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix1(2*k-3,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k) 
    matrix1(2*k-2,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix1(2*k-1,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix1(2*k,2*k-2) =   quad(@(x) d2_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix1(2*k+1,2*k-2) = quad(@(x) d2_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix1(2*k-4,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix1(2*k-3,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix1(2*k-2,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix1(2*k-1,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix1(2*k,2*k-1) =   quad(@(x) d2_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix1(2*k+1,2*k-1) = quad(@(x) d2_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis
matrix1(1,1) = quad(@(x) d2_phi_right{1}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)), startPoint, intervalLength+startPoint,Error)+alphan*quad(@(x) d2_psi_right{1}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error);
matrix1(2,1) = quad(@(x) d2_phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error)                            +alphan*quad(@(x) d2_psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error);
matrix1(3,1) = quad(@(x) d2_phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error)                            +alphan*quad(@(x) d2_psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error);


%Handle last interval 'half case'
%Phis 
matrix1(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) d2_phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betan*quad(@(x) d2_psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix1(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) d2_phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betan*quad(@(x) d2_psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix1(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) d2_phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(@(x) d2_psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);


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
matrix1(2*k-3,2*k-2) = quad(@(x) d2_phi_left{k}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) =   quad(@(x) d2_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix1(2*k+1,2*k-2) = quad(@(x) d2_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-3,2*k-1) = quad(@(x) d2_psi_left{k}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) = quad(@(x)   d2_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix1(2*k+1,2*k-1) = quad(@(x) d2_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix1(2*k-4,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix1(2*k-3,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) =   quad(@(x) d2_phi_right{k}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-4,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix1(2*k-3,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) =   quad(@(x) d2_psi_right{k}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

    Matrix1Finaln = Cn1*matrix1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of n''(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates for RHS of equation for p-type carriers (forcing function times test functions)
%This calculation is performed outside the iterative loop, as it will not
%change over the course of changing p(x), n(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RHS = Gopp;

%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(2*(numIntervals+1)-2,1);

for k = 2:numIntervals
    %Phi test function
    rhs_matrix(2*k-2) = quad(@(x) RHS(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) RHS(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi test function
    rhs_matrix(2*k-1) = quad(@(x) RHS(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) RHS(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    
end

%Handle 1st interval 'half case'
%Phi test function
rhs_matrix(1) = quad(@(x) RHS(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)),startPoint,intervalLength+startPoint, Error);    

%Handle last interval 'half case'
%Phi test function
rhs_matrix(2*(numIntervals+1)-2) = quad(@(x) RHS(x).*(phi_left{numIntervals+1}(x)+betap.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint,(numIntervals)*intervalLength+startPoint,Error);
   
RHSOutp =rhs_matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End RHS calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with p''(x) term
% %Performed outside iterative loop because this term is not influenced by
% %p(x) or n(x).  Should not change over calculations.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Front end coefficient of p''(x) term
 Cp1=ksi*T*mup/q - Dp;

%Set up empty matrix for holding coefficients
%Put Phi basis functions on odd rows, put Psi basis functions on even rows
%Subtract 2 rows and 2 columns to account for Boundary Condition
%limitations
 matrix1 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);
for k = 3:numIntervals-1
 
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix1(2*k-4,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix1(2*k-3,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k) 
    matrix1(2*k-2,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix1(2*k-1,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix1(2*k,2*k-2) =   quad(@(x) d2_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix1(2*k+1,2*k-2) = quad(@(x) d2_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix1(2*k-4,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix1(2*k-3,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix1(2*k-2,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix1(2*k-1,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix1(2*k,2*k-1) =   quad(@(x) d2_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix1(2*k+1,2*k-1) = quad(@(x) d2_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis
matrix1(1,1) = quad(@(x) d2_phi_right{1}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)), startPoint, intervalLength+startPoint,Error)+alphap*quad(@(x) d2_psi_right{1}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error);
matrix1(2,1) = quad(@(x) d2_phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error)                            +alphap*quad(@(x) d2_psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error);
matrix1(3,1) = quad(@(x) d2_phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error)                            +alphap*quad(@(x) d2_psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error);


%Handle last interval 'half case'
%Phis 
matrix1(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) d2_phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betap*quad(@(x) d2_psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix1(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) d2_phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betap*quad(@(x) d2_psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix1(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) d2_phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betap.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betap*quad(@(x) d2_psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betap.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);


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
matrix1(2*k-3,2*k-2) = quad(@(x) d2_phi_left{k}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are
%shifted up one spot. 

%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) =   quad(@(x) d2_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix1(2*k+1,2*k-2) = quad(@(x) d2_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-3,2*k-1) = quad(@(x) d2_psi_left{k}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) =   quad(@(x) d2_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix1(2*k+1,2*k-1) = quad(@(x) d2_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    
%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix1(2*k-4,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix1(2*k-3,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix1(2*k-2,2*k-2) = quad(@(x) d2_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix1(2*k-1,2*k-2) = quad(@(x) d2_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix1(2*k,2*k-2) =   quad(@(x) d2_phi_right{k}(x).*(phi_left{numIntervals+1}(x)+betap.*psi_left{numIntervals+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix1(2*k-4,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix1(2*k-3,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix1(2*k-2,2*k-1) = quad(@(x) d2_psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix1(2*k-1,2*k-1) = quad(@(x) d2_psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) d2_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix1(2*k,2*k-1) =   quad(@(x) d2_psi_right{k}(x).*(phi_left{numIntervals+1}(x)+betap.*psi_left{numIntervals+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

   Matrix1Finalp = Cp1*matrix1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of p''(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin iterative loop to calculate p(x) and n(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iterationCount=1:maxIterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate current representations of p(x),p'(x),p''(x)
%Values associated with 1st node in system are undefined
%Value defaults to 0, but is never used in code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_dist = cell(1,numIntervals+1);
pprime_dist =cell(1,numIntervals+1);
pdouble_dist =cell(1,numIntervals+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finds functions in the interval to the left of the kth node
%Note: This implies that p(1), p'(1), p''(1) are undefined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=2:numIntervals+1
p_dist{index}= matlabFunction(px(index,Cpjs,Dpjs,phi_right_syms,phi_left_syms,psi_right_syms,psi_left_syms));
pprime_dist{index}=matlabFunction(diff(px(index,Cpjs,Dpjs,phi_right_syms,phi_left_syms,psi_right_syms,psi_left_syms),x,1));
pdouble_dist{index}=matlabFunction(diff(px(index,Cpjs,Dpjs,phi_right_syms,phi_left_syms,psi_right_syms,psi_left_syms),x,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with n'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Front end coefficient for this matrix.  DOES NOT INCLUDE CONTRIBUTION FROM
%p(x) DISTRIBUTION.
     Cn2 = ksi*T*mun/q;

     matrix2 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%pterm = (p'(x)/p(x))
%term associated with phi_left(k) interval of n(x)
%pprime_dist{k}(x)./p_dist{k}(x)
%term associated with phi_right(k) interval of n(x)
%pprime_dist{k+1}(x)./p_dist{k+1}(x)
for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Cpk

    %Phi(k) Phi(k-1)    
    matrix2(2*k-4,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix2(2*k-3,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k)
    matrix2(2*k-2,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix2(2*k-1,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix2(2*k,2*k-2) =   quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix2(2*k+1,2*k-2) = quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix2(2*k-4,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix2(2*k-3,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix2(2*k-2,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix2(2*k-1,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix2(2*k,2*k-1) =   quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix2(2*k+1,2*k-1) = quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis 
matrix2(1,1) = quad(@(x) pprime_dist{2}(x)./p_dist{2}(x).*d_phi_right{1}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error)+alphan*quad(@(x) pprime_dist{2}(x)./p_dist{2}(x).*d_psi_right{1}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error);
matrix2(2,1) = quad(@(x) pprime_dist{2}(x)./p_dist{2}(x).*d_phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphan*quad(@(x) pprime_dist{2}(x)./p_dist{2}(x).*d_psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error);
matrix2(3,1) = quad(@(x) pprime_dist{2}(x)./p_dist{2}(x).*d_phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphan*quad(@(x) pprime_dist{2}(x)./p_dist{2}(x).*d_psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error);


%Handle last interval 'half case'
%Phis 
matrix2(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x).*d_phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betan*quad(@(x) pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x).*d_psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix2(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x).*d_phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betan*quad(@(x) pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x).*d_psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix2(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x).*d_phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan*quad(@(x) pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x).*d_psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);


%Handle interval 2 and n-1 case
%These intervals generate terms in the rows that would normally be in the
%for-loop matrix, but are omitted due to boundary conditions.  Due to the
%unusual number fo rows included in the final matrix, they must be input
%separately.

%Node 2
%Phis
%Code is copied for code in for-loop above
k=2;
%Phi(k) Phi(k-1)
matrix2(2*k-3,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Phi(k) Phi(k)
matrix2(2*k-2,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*phi_left{k}(x),       (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix2(2*k-1,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*psi_left{k}(x),       (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix2(2*k,2*k-2) =   quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix2(2*k+1,2*k-2) = quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix2(2*k-3,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix2(2*k-2,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*phi_left{k}(x),       (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix2(2*k-1,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*psi_left{k}(x),       (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix2(2*k,2*k-1) =   quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix2(2*k+1,2*k-1) = quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   

%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix2(2*k-4,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix2(2*k-3,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix2(2*k-2,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix2(2*k-1,2*k-2) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix2(2*k,2*k-2) =   quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_phi_right{k}(x).*(phi_left{k+1}(x)+betan.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix2(2*k-4,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix2(2*k-3,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix2(2*k-2,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix2(2*k-1,2*k-1) = quad(@(x) pprime_dist{k}(x)./p_dist{k}(x).*d_psi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix2(2*k,2*k-1) =   quad(@(x) pprime_dist{k+1}(x)./p_dist{k+1}(x).*d_psi_right{k}(x).*(phi_left{k+1}(x)+betan.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);

  Matrix2Finaln = Cn2*matrix2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of n'(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with n(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
matrix3 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%Coefficient Cn3 is built into filling the matrix since it has a difference
%that can't be removed built-in.

%%pterm = ksi*T*mun/q*(p''(x)/p(x)-(p'(x)/p(x))^2)-1/taun
%%term associated with phi_left(k) interval of n(x)
%%@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun)
%%term associated with phi_right(k) interval of n(x)
%%@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun)

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix3(2*k-4,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix3(2*k-3,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k)
    matrix3(2*k-2,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix3(2*k-1,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix3(2*k,2*k-2) =   quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix3(2*k+1,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix3(2*k-4,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix3(2*k-3,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix3(2*k-2,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix3(2*k-1,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix3(2*k,2*k-1) =   quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix3(2*k+1,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis
matrix3(1,1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{2}(x)./p_dist{2}(x) - (pprime_dist{2}(x)./p_dist{2}(x)).^2)-1/taun).*phi_right{1}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error)+alphan*quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{2}(x)./p_dist{2}(x) - (pprime_dist{2}(x)./p_dist{2}(x)).^2)-1/taun).*psi_right{1}(x).*(phi_right{1}(x)+alphan.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error);
matrix3(2,1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{2}(x)./p_dist{2}(x) - (pprime_dist{2}(x)./p_dist{2}(x)).^2)-1/taun).*phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphan*quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{2}(x)./p_dist{2}(x) - (pprime_dist{2}(x)./p_dist{2}(x)).^2)-1/taun).*psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error);
matrix3(3,1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{2}(x)./p_dist{2}(x) - (pprime_dist{2}(x)./p_dist{2}(x)).^2)-1/taun).*phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphan*quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{2}(x)./p_dist{2}(x) - (pprime_dist{2}(x)./p_dist{2}(x)).^2)-1/taun).*psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error);

%Handle last interval 'half case'
%Phis
matrix3(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x) - (pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x)).^2)-1/taun).*phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betan.*quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x) - (pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x)).^2)-1/taun).*psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix3(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x) - (pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x)).^2)-1/taun).*phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betan.*quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x) - (pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x)).^2)-1/taun).*psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix3(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x) - (pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x)).^2)-1/taun).*phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betan.*quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x) - (pprime_dist{numIntervals+1}(x)./p_dist{numIntervals+1}(x)).^2)-1/taun).*psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betan.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);

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
matrix3(2*k-3,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*(phi_right{k-1}(x)+alphan.*psi_right{k-1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Phi(k) Phi(k)
matrix3(2*k-2,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*phi_left{k}(x),           (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix3(2*k-1,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*psi_left{k}(x),           (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
 matrix3(2*k,2*k-2)   = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
 matrix3(2*k+1,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix3(2*k-3,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*(phi_right{k-1}(x)+alphan.*psi_right{k-1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix3(2*k-2,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*phi_left{k}(x),           (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix3(2*k-1,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*psi_left{k}(x),           (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
 matrix3(2*k,2*k-1) =   quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
 matrix3(2*k+1,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix3(2*k-4,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
matrix3(2*k-3,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix3(2*k-2,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix3(2*k-1,2*k-2) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*phi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix3(2*k,2*k-2) =   quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*phi_right{k}(x).*(phi_left{k+1}(x)+betan.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix3(2*k-4,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix3(2*k-3,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix3(2*k-2,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix3(2*k-1,2*k-1) = quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k}(x)./p_dist{k}(x) - (pprime_dist{k}(x)./p_dist{k}(x)).^2)-1/taun).*psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix3(2*k,2*k-1) =   quad(@(x) (ksi.*T.*mun./q.*(pdouble_dist{k+1}(x)./p_dist{k+1}(x) - (pprime_dist{k+1}(x)./p_dist{k+1}(x)).^2)-1/taun).*psi_right{k}(x).*(phi_left{k+1}(x)+betan.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);


   Matrix3Finaln = matrix3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of n(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%Combine expressions from each portion for total matrix of LHS
LHSn = Matrix1Finaln +Matrix2Finaln +Matrix3Finaln;



%Solve simultaneous equations for basis function coefficients
Coeffsn = LHSn\RHSOutn;

%Collect Cj's in an array for display purposes
Cnjs = zeros(1,numIntervals+1);
Dnjs = zeros(1,numIntervals+1);

for index=2:length(Cpjs)-1
   Cnjs(1,index) = Coeffsn(2*index-2);
   Dnjs(1,index) = Coeffsn(2*index-1);
end

Cnjs(1,1) = Coeffsn(1);
Dnjs(1,1) = alphan*Coeffsn(1);
Cnjs(1,numIntervals+1) = Coeffsn(length(Coeffsn));
Dnjs(1,numIntervals+1) = betan*Coeffsn(length(Coeffsn));


fig4 = figure;
plot(indexn,Cnjs,indexn,subs(nGuess,x,indexn));title('n(x) calculated');
saveas(fig4,strcat(filepath,'/n_calculated_iteration',num2str(iterationCount),'.jpg'))

fig6 = figure;
plot(indexn,Dnjs,indexn,subs(diff(nGuess,x,1),x,indexn)),title('nprime calculated');
saveas(fig6,strcat(filepath,'/nprime_calculated_iteration',num2str(iterationCount),'.jpg'))

%Collect Cj's in an array for display purposes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up n(x), n'(x), and n''(x) distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_dist = cell(1,numIntervals+1);
nprime_dist =cell(1,numIntervals+1);
ndouble_dist =cell(1,numIntervals+1);

for index=2:numIntervals+1
n_dist{index}= matlabFunction(nx(index,Cnjs,Dnjs,phi_right_syms,phi_left_syms,psi_right_syms,psi_left_syms));
nprime_dist{index}=matlabFunction(diff(nx(index,Cnjs,Dnjs,phi_right_syms,phi_left_syms,psi_right_syms,psi_left_syms),x,1));
ndouble_dist{index}=matlabFunction(diff(nx(index,Cnjs,Dnjs,phi_right_syms,phi_left_syms,psi_right_syms,psi_left_syms),x,2));
end


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Begin calculation of p-type carriers
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with p'(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Front end coefficient for this matrix.  
Cp2 = -ksi*T*mup/q;
 
matrix2 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%%pterm = (p'(x)/p(x))
%%term associated with phi_left(k) interval of n(x)
%%nprime_dist{k}(x)./n_dist{k}(x)
%%term associated with phi_right(k) interval of n(x)
%%nprime_dist{k+1}(x)./n_dist{k+1}(x)

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Cpk
    
    %Phi(k) Phi(k-1)    
    matrix2(2*k-4,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix2(2*k-3,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k)
    matrix2(2*k-2,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix2(2*k-1,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix2(2*k,2*k-2) =   quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix2(2*k+1,2*k-2) = quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix2(2*k-4,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix2(2*k-3,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix2(2*k-2,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*phi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix2(2*k-1,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*psi_left{k}(x),    (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix2(2*k,2*k-1) =   quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix2(2*k+1,2*k-1) = quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis
matrix2(1,1) = quad(@(x) nprime_dist{2}(x)./n_dist{2}(x).*d_phi_right{1}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error)+alphap*quad(@(x) nprime_dist{2}(x)./n_dist{2}(x).*d_psi_right{1}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error);
matrix2(2,1) = quad(@(x) nprime_dist{2}(x)./n_dist{2}(x).*d_phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphap*quad(@(x) nprime_dist{2}(x)./n_dist{2}(x).*d_psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error);
matrix2(3,1) = quad(@(x) nprime_dist{2}(x)./n_dist{2}(x).*d_phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphap*quad(@(x) nprime_dist{2}(x)./n_dist{2}(x).*d_psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error);


%Handle last interval 'half case'
%Phis
matrix2(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x).*d_phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betap*quad(@(x) nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x).*d_psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix2(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x).*d_phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                      +betap*quad(@(x) nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x).*d_psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix2(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x).*d_phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betap.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betap*quad(@(x) nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x).*d_psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+beta.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);

%Handle interval 2 and n-1 case
%These intervals generate terms in the rows that would normally be in the
%for-loop matrix, but are omitted due to boundary conditions.  Due to the
%unusual number fo rows included in the final matrix, they must be input
%separately.

%Node 2
%Phis
%Code is copied for code in for-loop above
k=2;
%Phi(k) Phi(k-1)
matrix2(2*k-3,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*(phi_right{k-1}(x)+alphap.*psi_right{k-1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Phi(k) Phi(k)
matrix2(2*k-2,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*phi_left{k}(x),       (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix2(2*k-1,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*psi_left{k}(x),       (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix2(2*k,2*k-2) =   quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix2(2*k+1,2*k-2) = quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix2(2*k-3,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*(phi_right{k-1}(x)+alphap.*psi_right{k-1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix2(2*k-2,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*phi_left{k}(x),       (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix2(2*k-1,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*psi_left{k}(x),       (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix2(2*k,2*k-1) =   quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix2(2*k+1,2*k-1) = quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   

%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix2(2*k-4,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
%Omitted because of Psi(k-1) test function term, all subsequent terms are shifted up one spot. 
matrix2(2*k-3,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix2(2*k-2,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix2(2*k-1,2*k-2) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_phi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix2(2*k,2*k-2) =   quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_phi_right{k}(x).*(phi_left{k+1}(x)+betap.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix2(2*k-4,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix2(2*k-3,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix2(2*k-2,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix2(2*k-1,2*k-1) = quad(@(x) nprime_dist{k}(x)./n_dist{k}(x).*d_psi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix2(2*k,2*k-1) =   quad(@(x) nprime_dist{k+1}(x)./n_dist{k+1}(x).*d_psi_right{k}(x).*(phi_left{k+1}(x)+betap*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);

  Matrix2Finalp = Cp2*matrix2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of p'(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with p(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
   matrix3 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%Coefficient Cp3 is built into filling the matrix since it has a difference
%that can't be removed built-in.

%%pterm = ksi*T*mup/q*((n'(x)/n(x))^2-n''(x)/n(x))-1/taup
%%term associated with phi_left(k) interval of p(x)
%%@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup)
%%term associated with phi_right(k) interval of p(x)
%%@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup)

for k = 3:numIntervals-1
   
    %%Calculate for interval's Phi basis function
    %%Corresponds to column Ck
    %Phi(k) Phi(k-1)
    matrix3(2*k-4,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Phi(k) Psi(k-1)
    matrix3(2*k-3,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Phi(k) Phi(k)
    matrix3(2*k-2,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k)
    matrix3(2*k-1,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Phi(k+1)
    matrix3(2*k,2*k-2) =   quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Phi(k) Psi(k+1)
    matrix3(2*k+1,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);   
     

    %%Calculate for interval's Psi basis function
    %Corresponds to column Dk
    %%Psi(k) Phi(k-1)
    matrix3(2*k-4,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

    %Psi(k) Psi(k-1)
    matrix3(2*k-3,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
    %Psi(k) Phi(k)
    matrix3(2*k-2,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k)
    matrix3(2*k-1,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
    %Psi(k) Phi(k+1)
    matrix3(2*k,2*k-1) =   quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
    %Psi(k) Psi(k+1)
    matrix3(2*k+1,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
end

%Handle 1st interval 'half case'
%Phis

matrix3(1,1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{2}(x)./n_dist{2}(x)).^2 - ndouble_dist{2}(x)./n_dist{2}(x))-1/taup).*phi_right{1}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error)+alphap*quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{2}(x)./n_dist{2}(x)).^2 - ndouble_dist{2}(x)./n_dist{2}(x))-1/taup).*psi_right{1}(x).*(phi_right{1}(x)+alphap.*psi_right{1}(x)),startPoint, intervalLength+startPoint,Error);
matrix3(2,1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{2}(x)./n_dist{2}(x)).^2 - ndouble_dist{2}(x)./n_dist{2}(x))-1/taup).*phi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphap*quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{2}(x)./n_dist{2}(x)).^2 - ndouble_dist{2}(x)./n_dist{2}(x))-1/taup).*psi_right{1}(x).*phi_left{2}(x),startPoint, intervalLength+startPoint,Error);
matrix3(3,1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{2}(x)./n_dist{2}(x)).^2 - ndouble_dist{2}(x)./n_dist{2}(x))-1/taup).*phi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error)                           +alphap*quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{2}(x)./n_dist{2}(x)).^2 - ndouble_dist{2}(x)./n_dist{2}(x))-1/taup).*psi_right{1}(x).*psi_left{2}(x),startPoint, intervalLength+startPoint,Error);

%Handle last interval 'half case'
%Phis

matrix3(2*(numIntervals)-2,2*(numIntervals+1)-2)= quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x)).^2 - ndouble_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x))-1/taup).*phi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                     +betap*quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x)).^2 - ndouble_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x))-1/taup).*psi_left{numIntervals+1}(x).*phi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix3(2*(numIntervals)-1,2*(numIntervals+1)-2)= quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x)).^2 - ndouble_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x))-1/taup).*phi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)                                     +betap*quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x)).^2 - ndouble_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x))-1/taup).*psi_left{numIntervals+1}(x).*psi_right{numIntervals}(x),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);
matrix3(2*(numIntervals),2*(numIntervals+1)-2)=   quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x)).^2 - ndouble_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x))-1/taup).*phi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betap*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error)+betap*quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x)).^2 - ndouble_dist{numIntervals+1}(x)./n_dist{numIntervals+1}(x))-1/taup).*psi_left{numIntervals+1}(x).*(phi_left{numIntervals+1}(x)+betap.*psi_left{numIntervals+1}(x)),(numIntervals-1)*intervalLength+startPoint, (numIntervals)*intervalLength+startPoint,Error);


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
matrix3(2*k-3,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*(phi_right{k-1}(x)+alphap.*psi_right{k-1}(x)),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Phi(k) Phi(k)
matrix3(2*k-2,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*phi_left{k}(x),          (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix3(2*k-1,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*psi_left{k}(x),          (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix3(2*k,2*k-2)   = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k+1)
matrix3(2*k+1,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);  


%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix3(2*k-3,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*(phi_right{k-1}(x)+alphap.*psi_right{k-1}(x)), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
  
%Psi(k) Phi(k)
matrix3(2*k-2,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*phi_left{k}(x),          (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix3(2*k-1,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*psi_left{k}(x),          (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix3(2*k,2*k-1) =   quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*phi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Psi(k) Psi(k+1)
matrix3(2*k+1,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*psi_left{k+1}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%Interval n-1
%Phis
%Code is copied for code in for-loop above
k=numIntervals;
%Phi(k) Phi(k-1)
matrix3(2*k-4,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*phi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Psi(k-1)
matrix3(2*k-3,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*psi_right{k-1}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
 
%Phi(k) Phi(k)
matrix3(2*k-2,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*phi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Psi(k)
matrix3(2*k-1,2*k-2) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*phi_left{k}(x).*psi_left{k}(x),   (k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    
%Phi(k) Phi(k+1)
matrix3(2*k,2*k-2) =   quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*phi_right{k}(x).*(phi_left{k+1}(x)+betap.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
    

%%Calculate for interval's Psi basis function
%Corresponds to column Dk
%%Psi(k) Phi(k-1)
matrix3(2*k-4,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*phi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);

%Psi(k) Psi(k-1)
matrix3(2*k-3,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*psi_right{k-1}(x), (k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error);
   
%Psi(k) Phi(k)
matrix3(2*k-2,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*phi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*phi_right{k}(x),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);
   
%Psi(k) Psi(k)
matrix3(2*k-1,2*k-1) = quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k}(x)./n_dist{k}(x)).^2 - ndouble_dist{k}(x)./n_dist{k}(x))-1/taup).*psi_left{k}(x).*psi_left{k}(x),(k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,Error)+quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*psi_right{k}(x),(k-1)*intervalLength+startPoint,(k)*intervalLength+startPoint,Error);
    
%Psi(k) Phi(k+1)
matrix3(2*k,2*k-1) =   quad(@(x) (ksi.*T.*mup./q.*((nprime_dist{k+1}(x)./n_dist{k+1}(x)).^2 -ndouble_dist{k+1}(x)./n_dist{k+1}(x))-1/taup).*psi_right{k}(x).*(phi_left{k+1}(x)+betap.*psi_left{k+1}(x)),(k-1)*intervalLength+startPoint, (k)*intervalLength+startPoint,Error);


Matrix3Finalp = matrix3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of p(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Combine expressions from each portion for total matrix of LHS
LHSp = Matrix1Finalp +Matrix2Finalp +Matrix3Finalp;



%Solve simultaneous equations for basis function coefficients
Coeffsp = LHSp\RHSOutp;

%Collect Cj's in an array for display purposes
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);


for index=2:length(Cpjs)-1
   Cpjs(1,index) = Coeffsp(2*index-2);
   Dpjs(1,index) = Coeffsp(2*index-1);
end

Cpjs(1,1) = Coeffsp(1);
Dpjs(1,1) = alphap*Coeffsp(1);

Cpjs(1,numIntervals+1) = Coeffsp(length(Coeffsp));
Dpjs(1,numIntervals+1) = betap*Coeffsp(length(Coeffsp));



fig5 = figure;
plot(indexp,Cpjs,indexp,subs(pGuess,x,indexp));title('p(x) calculated');
saveas(fig5,strcat(filepath,'/p_calculated_iteration iteration',num2str(iterationCount),'.jpg'))


fig7 = figure;
plot(indexp,Dpjs,indexp,subs(diff(pGuess,x,1),x,indexp)),title('pprime calculated');
saveas(fig7,strcat(filepath,'/pprime_calculated_iteration iteration',num2str(iterationCount),'.jpg'))


save(strcat(filepath,'/','iteration',num2str(iterationCount),'.mat'))


iterationCount

end

% profile viewer
% p=profile('info')
% profile off
% profsave(p,strcat(filepath,'/profile_results'))


%Displays runtime of m-file
time = cputime-t

%Store variable values to be used later
save(strcat(filepath,'/',num2str(numIntervals),'terms.mat'))