clear all
close all
t=cputime;
%profile on

%Modified Newton, noisy input, Jacobian, terms one only


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
ksiN = 1.38*10^-23;
ksiP = 1.38*10^-23;
T = 300; %(*Cell temperature in K*)
h = 6.626 * 10^-34;
Ec = -chi;
Ev = -gamma;
e0m = 8.854*10^-15;  %(*permittivity of free space in [s^4 A^2/m^3 kg]*)
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
dx = intervalLength;
maxIterations=5;
% Store node locations
node = startPoint:dx:totalLength;

%Random noise bounds
rand=0.8;

AbsError=10^-18;
RelError=10^-6;

filepath=strcat('modifiedNewton v6/rand',num2str(rand),'/',num2str(numIntervals),'terms');
mkdir(filepath)

%(*Solar incidence data*)
alphaOpt = 10000; %(*/cm*) (*GaAs*)
Conc = 1;
Lflux = 2.0*10^17; %(*/cm^2 s *) (*photons absorbed by GaAs at 1 sun*)
R = 0.1; %(*reflection coefficient*)
Gop = @(x) alphaOpt*Conc*Lflux*(1 - R)*exp(-alphaOpt*(x));%(*/cm^3*) 


indexp = linspace(startPoint,totalLength+startPoint,numIntervals+1);
indexn = linspace(startPoint,totalLength+startPoint,numIntervals+1);

%Define general basis functions for interval 0 to 1
% Use inline functions; define generic basis functions and derivatives

%Phi function is chosen such that phi_right (f1) has f(0) = 1, f(1)=0, f'(0) = 1,
%f'(1)=0
%Phi function is chosen such that phi_left (f2) has f(0) = 0, f(1)=1, f'(0) = 1,
%f'(1)=0
%f1 = phi_right
%f2 = phi_left
f1 =@(x) (1-3.*(x/dx).^2+ 2.*(x/dx).^3);
f2 =@(x) (3*(x/dx).^2-2*(x/dx).^3); 

f1prime = @(x) (-6/dx^2.*(x) + 6/dx^3.*(x).^2);
f2prime = @(x) (6/dx^2.*(x) - 6/dx^3.*(x).^2);

f1primeprime = @(x) (-6/dx^2 + 12/dx^3.*(x));
f2primeprime = @(x) (6/dx^2 - 12/dx^3.*(x));

%Psi function is chosen such that psi_right (f3) has f'(0) = 1, f'(1)=0, f(0) =
%0, f(1)=0
%Psi function is chosen such that psi_left (f4) has f'(0) = 0, f'(1)=1, f(0) = 0,
%f(1)=0
%f3 = psi_left
%f4 = psi_right
f3 = @(x) ((x/dx)- 2*(x/dx).^2+(x/dx).^3)*dx;
f4 = @(x) (-(x/dx).^2+(x/dx).^3)*dx;


f3prime = @(x) (1/dx - 4/dx.^2.*(x) + 3/dx.^3.*(x).^2)*dx;
f4prime = @(x) (-2/dx^2.*(x) + 3/dx^3.*(x).^2)*dx;

f3primeprime = @(x) (-4/dx.^2 + 6/dx.^3.*(x))*dx;
f4primeprime = @(x) (-2/dx.^2 + 6/dx.^3.*(x))*dx;

% x_axis=0:dx/100:dx;
% subplot(221),plot(x_axis,f1(x_axis))
% subplot(222),plot(x_axis,f2(x_axis))
% subplot(223),plot(x_axis,f3(x_axis))
% subplot(224),plot(x_axis,f4(x_axis))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up known n(x) and p(x) distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Set up initial guess for p(x).  The curve goes has the vague shape of the
%expected output (peak near 0 and vaguely exponential).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: If pGuess is change, initial BCs for pGuess must also be changed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1=8.14812*10^13;
p2=-2.19257*10^11;
p3=-8.12149*10^13;
p4=1.33755*10^12;
gammap1=100;
gammap2=100000;
gammap3=1;
gammap4=10;

%Symbolic form of pGuess
pGuessOrig = @(x) p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(-gammap3*x)+p4*exp(-gammap4*x);

n1=-1.79768*10^10;
n2=1.29154*10^11;
n3=9.94351*10^11;
n4=5.6616*10^9;
gamman1=100000;
gamman2=10000;
gamman3=1;
gamman4=10;

%Symbolic form of nGuess
nGuessOrig = @(x) n1*exp(-gamman1*x)+n2*exp(-gamman2*x)+n3*exp(-gamman3*x)+n4*exp(-gamman4*x);

nprimeOrig = @(x) -(n1*gamman1*exp(-gamman1*x) + n2*gamman2*exp(-gamman2*x) + n3*gamman3*exp(-gamman3*x)+n4*gamman4*exp(-gamman4*x));
pprimeOrig = @(x) -(p1*gammap1*exp(-gammap1*x) + p2*gammap2*exp(-gammap2*x) + p3*gammap3*exp(-gammap3*x)+p4*gammap4*exp(-gammap4*x));

nprimeprimeOrig = @(x) (n1*gamman1*gamman1*exp(-gamman1*x) + n2*gamman2*gamman2*exp(-gamman2*x) + n3*gamman3*gamman3*exp(-gamman3*x)+n4*gamman4*gamman4*exp(-gamman4*x));
pprimeprimeOrig = @(x) (p1*gammap1*gammap1*exp(-gammap1*x) + p2*gammap2*gammap2*exp(-gammap2*x) + p3*gammap3*gammap3*exp(-gammap3*x)+p4*gammap4*gammap4*exp(-gammap4*x));

cn1 = Dn - ksiN*T*mun/q;
cn2 = ksiN*T*mun/q;

cp1 = ksiP*T*mup/q - Dp;
cp2 = -ksiP*T*mup/q;

%Check BCs
% alphaNtest = nprimeOrig(0)/nGuessOrig(0)
% alphaPtest = pprimeOrig(0)/pGuessOrig(0)
% betaNtest = nprimeOrig(Wn)/nGuessOrig(Wn)
% betaPtest = pprimeOrig(Wn)/pGuessOrig(Wn)


%Gopn = @(x) cn1.*(n1*(gamman1*gamman1)*exp(-gamman1*x) + n2*(gamman2*gamman2)*exp(-gamman2*x) + n3*(gamman3*gamman3)*exp(-gamman3*x))    +(ksiN*T*mun/q*(pprimeprimeOrig(x)./pGuessOrig(x) - (pprimeOrig(x)./pGuessOrig(x)).^2) - 1/taun).*(n1*exp(-gamman1*x) + n2*exp(-gamman2*x) + n3*exp(-gamman3*x));%+cn2.*pprime(x)./pGuess(x).*(-n1*(gamman1)*exp(-gamman1*x) - n2*(gamman2)*exp(-gamman2*x) - n3*(gamman3)*exp(-gamman3*x))     
%Gopp = @(x) cp1.*(p1*(gammap1*gammap1)*exp(-gammap1*x) + p2*(gammap2*gammap2)*exp(-gammap2*x) + p3*(gammap3*gammap3)*exp(-gammap3*x))    +(ksiP*T*mup/q*((nprimeOrig(x)./nGuessOrig(x)).^2 - nprimeprimeOrig(x)./nGuessOrig(x)) - 1/taup).*(p1*exp(-gammap1*x) + p2*exp(-gammap2*x) + p3*exp(-gammap3*x));%+cp2.*nprime(x)./nGuess(x).*(-p1*(gammap1)*exp(-gammap1*x) - p2*(gammap2)*exp(-gammap2*x) - p3*(gammap3)*exp(-gammap3*x))     

%Adjust pGuess to include error
%Adjust pGuess to include error

%'Noisified' terms
%n3*rand,gamman1,gamman2,gamman3 set;n1 and n2 calculated from those
%numbers using BCs
%Calcs done in Mathematica

n1rand=-2.43711*10^10;
n2rand=1.53127*10^11;
n3rand=1.07209*10^13;
n4rand=-9.73843*10^12;
gamman1rand=100000*rand;

nGuess = @(x) (n1rand*exp(-gamman1rand*x)+n2rand*exp(-gamman2*x)+n3rand*exp(-gamman3*x)+n4rand*exp(-gamman4*x));
nprime = @(x) -(n1rand*gamman1rand*exp(-gamman1rand*x) + n2rand*gamman2*exp(-gamman2*x) + n3rand*gamman3*exp(-gamman3*x)+n4rand*gamman4*exp(-gamman4*x));
nprimeprime = @(x) n1rand*gamman1rand*gamman1rand*exp(-gamman1rand*x) + n2rand*gamman2*gamman2*exp(-gamman2*x) + n3rand*gamman3*gamman3*exp(-gamman3*x)+n4rand*gamman4*gamman4*exp(-gamman4*x);

%p3*rand,gammap1,gammap2,gammap3 set;p1 and p2 calculated from those
%numbers using BCs
%Calcs done in Mathematica

p1rand=1.71267*10^15;
p2rand=-2.92243*10^11;
p3rand=1.60693*10^16;
p4rand=-1.77803*10^16;
gammap2rand=100000*rand;

pGuess = @(x) (p1rand*exp(-gammap1*x)+p2rand*exp(-gammap2rand*x)+p3rand*exp(-gammap3*x)+p4rand*exp(-gammap4*x));
pprime = @(x) -(p1rand*gammap1*exp(-gammap1*x) + p2rand*gammap2rand*exp(-gammap2rand*x) + p3rand*gammap3*exp(-gammap3*x)+p4rand*gammap4*exp(-gammap4*x));
pprimeprime =@(x) p1rand*gammap1*gammap1*exp(-gammap1*x) + p2rand*gammap2rand*gammap2rand*exp(-gammap2rand*x) + p3rand*gammap3*gammap3*exp(-gammap3*x)+p4rand*gammap4*gammap4*exp(-gammap4*x);

% alphaNtestNoise = nprime(0)/nGuess(0)
% alphaPtestNoise = pprime(0)/pGuess(0)
% betaNtestNoise = nprime(Wn)/nGuess(Wn)
% betaPtestNoise = pprime(Wn)/pGuess(Wn)
% 
% diffAlphaN=alphaNtest-alphaNtestNoise
% diffAlphaP=alphaPtest-alphaPtestNoise
% diffBetaN=betaNtest-betaNtestNoise
% diffBetaP=betaPtest-betaPtestNoise

% fig2 = figure;
% plot(node,nGuess(node)), title('nGuess')
% saveas(fig2,strcat(filepath,'/nGuess.jpg'))

%Convert initial guess into form finite element composition
%This will allow it to be usable in later calculations

%Boundary conditions of test function put into form usable for the
%calculations below
 
temp1 = ksiN*T*mun/q*(pprimeprime(node)./pGuess(node) - (pprime(node)./pGuess(node)).^2) -1/taun;
temp2 = ksiN*T*mun/q*(pprimeprimeOrig(node)./pGuessOrig(node) - (pprimeOrig(node)./pGuessOrig(node)).^2)-1/taun;

plot(node,temp1,node,temp2,'x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate composition of n(x) in finite element space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alpha is boundary conditions term, n'(0)=alpha*n(0)
alphan = nprimeOrig(0)/nGuessOrig(0);

%beta is boundary condition term, n'(L) = beta n(L)
betan= nprimeOrig(totalLength)/nGuessOrig(totalLength);

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alphap = pprimeOrig(0)/pGuessOrig(0);

%beta is boundary condition term, n'(L) = beta n(L)
betap= pprimeOrig(totalLength)/pGuessOrig(totalLength);


massmatrix = sparse (4*numIntervals,4*numIntervals);

massrhs = zeros(4*numIntervals,1);
for k = 1:numIntervals-2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 1  - (Upper left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %massmatrix(phi(n,r),xx)
	massmatrix(4*k-1,4*k-1) = massmatrix(4*k-1,4*k-1) + integral(@(x) f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k,4*k-1)   = massmatrix(4*k,4*k-1) +   integral(@(x) f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(n,r),xx)
    massmatrix(4*k-1,4*k) = massmatrix(4*k-1,4*k) + integral(@(x) f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k,4*k)   = massmatrix(4*k,4*k) +   integral(@(x) f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(phi(p,r),xx)
	massmatrix(4*k+1,4*k+1) = massmatrix(4*k+1,4*k+1) + integral(@(x) f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+2,4*k+1) = massmatrix(4*k+2,4*k+1) + integral(@(x) f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(p,r),xx)
	massmatrix(4*k+1,4*k+2) = massmatrix(4*k+1,4*k+2) + integral(@(x) f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+2,4*k+2) = massmatrix(4*k+2,4*k+2) + integral(@(x) f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2  - (Upper right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %massmatrix(phi(n,r),xx)
	massmatrix(4*k-1,4*k+3) = massmatrix(4*k-1,4*k+3) + integral(@(x) f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k,4*k+3)   = massmatrix(4*k,4*k+3) +   integral(@(x) f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(n,r),xx)
    massmatrix(4*k-1,4*k+4) = massmatrix(4*k-1,4*k+4) + integral(@(x) f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k,4*k+4)   = massmatrix(4*k,4*k+4) +   integral(@(x) f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(phi(p,r),xx)
	massmatrix(4*k+1,4*k+5) = massmatrix(4*k+1,4*k+5) + integral(@(x) f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+2,4*k+5) = massmatrix(4*k+2,4*k+5) + integral(@(x) f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(p,r),xx)
	massmatrix(4*k+1,4*k+6) = massmatrix(4*k+1,4*k+6) + integral(@(x) f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+2,4*k+6) = massmatrix(4*k+2,4*k+6) + integral(@(x) f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 3  - (Bottom left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %massmatrix(phi(n,r),xx)
	massmatrix(4*k+3,4*k-1) = massmatrix(4*k+3,4*k-1) + integral(@(x) f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+4,4*k-1) = massmatrix(4*k+4,4*k-1) + integral(@(x) f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(n,r),xx)
    massmatrix(4*k+3,4*k) = massmatrix(4*k+3,4*k) + integral(@(x) f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+4,4*k) = massmatrix(4*k+4,4*k) + integral(@(x) f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(phi(p,r),xx)
	massmatrix(4*k+5,4*k+1) = massmatrix(4*k+5,4*k+1) + integral(@(x) f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+6,4*k+1) = massmatrix(4*k+6,4*k+1) + integral(@(x) f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(p,r),xx)
	massmatrix(4*k+5,4*k+2) = massmatrix(4*k+5,4*k+2) + integral(@(x) f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+6,4*k+2) = massmatrix(4*k+6,4*k+2) + integral(@(x) f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 4  - (Bottom right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %massmatrix(phi(n,r),xx)
	massmatrix(4*k+3,4*k+3) = massmatrix(4*k+3,4*k+3) + integral(@(x) f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+4,4*k+3) = massmatrix(4*k+4,4*k+3) + integral(@(x) f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(n,r),xx)
    massmatrix(4*k+3,4*k+4) = massmatrix(4*k+3,4*k+4) + integral(@(x) f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+4,4*k+4) = massmatrix(4*k+4,4*k+4) + integral(@(x) f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(phi(p,r),xx)
	massmatrix(4*k+5,4*k+5) = massmatrix(4*k+5,4*k+5) + integral(@(x) f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+6,4*k+5) = massmatrix(4*k+6,4*k+5) + integral(@(x) f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %massmatrix(psi(p,r),xx)
	massmatrix(4*k+5,4*k+6) = massmatrix(4*k+5,4*k+6) + integral(@(x) f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massmatrix(4*k+6,4*k+6) = massmatrix(4*k+6,4*k+6) + integral(@(x) f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RHS input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %nGuess terms
	massrhs(4*k-1) =      massrhs(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massrhs(4*k) =        massrhs(4*k)   +        integral(@(x) f3((x-node(k+1))).*nGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    massrhs(4*k+3) =      massrhs(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massrhs(4*k+4) =      massrhs(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %pGuess terms
	massrhs(4*k+1) =    massrhs(4*k+1) +      integral(@(x) f1((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massrhs(4*k+2) =    massrhs(4*k+2) +      integral(@(x) f3((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    massrhs(4*k+5) =    massrhs(4*k+5) +      integral(@(x) f2((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	massrhs(4*k+6) =    massrhs(4*k+6) +      integral(@(x) f4((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
end    
    
%Handle 1st interval 'half case'
%Top left 2x2 block
massmatrix(1,1)= integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrix(2,2)= integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Left 2x4 block
massmatrix(3,1) = integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrix(4,1) = integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

massmatrix(5,2) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrix(6,2) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Right 4x2 block
massmatrix(1,3) = integral(@(x) f2((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrix(1,4) = integral(@(x) f4((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

massmatrix(2,5) = integral(@(x) f2((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrix(2,6) = integral(@(x) f4((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%massmatrix(phi(n,r),xx)
massmatrix(3,3) = massmatrix(3,3) + integral(@(x) f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(4,3)   = massmatrix(4,3) + integral(@(x) f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

%massmatrix(psi(n,r),xx)
massmatrix(3,4) = massmatrix(3,4) + integral(@(x) f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(4,4)   = massmatrix(4,4) + integral(@(x) f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

%massmatrix(phi(p,r),xx)
massmatrix(5,5) = massmatrix(5,5) + integral(@(x) f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(6,5) = massmatrix(6,5) + integral(@(x) f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

%massmatrix(psi(p,r),xx)
massmatrix(5,6) = massmatrix(5,6) + integral(@(x) f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(6,6) = massmatrix(6,6) + integral(@(x) f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

%RHS nGuess edge cases
massrhs(1) = massrhs(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*nGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massrhs(3) = massrhs(3) + integral(@(x) f2((x-node(1))).*nGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massrhs(4) = massrhs(4) + integral(@(x) f4((x-node(1))).*nGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%RHS nGuess edge cases
massrhs(2) = massrhs(2) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*pGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massrhs(5) = massrhs(5) + integral(@(x) f2((x-node(1))).*pGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massrhs(6) = massrhs(6) + integral(@(x) f4((x-node(1))).*pGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Handle last interval 'half case'
last = 4*numIntervals;

%Handle 1st interval 'half case'
%Bottom right 2x2 block
massmatrix(last,last) =     integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-1,last-1) = integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Bottom 2x4 block
massmatrix(last,last-2) = integral(@(x) f3((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last,last-3) = integral(@(x) f1((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

massmatrix(last-1,last-4) = integral(@(x) f3((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-1,last-5) = integral(@(x) f1((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Left 4x2 Block
massmatrix(last-2,last) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-3,last) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

massmatrix(last-4,last-1) = integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-5,last-1) = integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%massmatrix(phi(n,r),xx)
massmatrix(last-5,last-5) = massmatrix(last-5,last-5) + integral(@(x) f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-4,last-5) = massmatrix(last-4,last-5) + integral(@(x) f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

%massmatrix(psi(n,r),xx)
massmatrix(last-5,last-4) = massmatrix(last-5,last-4) + integral(@(x) f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-4,last-4) = massmatrix(last-4,last-4) + integral(@(x) f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

%massmatrix(phi(p,r),xx)
massmatrix(last-3,last-3) = massmatrix(last-3,last-3) + integral(@(x) f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-2,last-3) = massmatrix(last-2,last-3) + integral(@(x) f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

%massmatrix(psi(p,r),xx)
massmatrix(last-3,last-2) = massmatrix(last-3,last-2) + integral(@(x) f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
massmatrix(last-2,last-2) = massmatrix(last-2,last-2) + integral(@(x) f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
%RHS edge cases
%nGuess
massrhs(last-1) = massrhs(last)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*nGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massrhs(last-5) = massrhs(last-5) + integral(@(x) f1((x-node(end)+dx)).*nGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massrhs(last-4) = massrhs(last-4) + integral(@(x) f3((x-node(end)+dx)).*nGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%pGuess
massrhs(last) = massrhs(last)+ integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*pGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massrhs(last-3) = massrhs(last-3) + integral(@(x) f1((x-node(end)+dx)).*pGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massrhs(last-2) = massrhs(last-2) + integral(@(x) f3((x-node(end)+dx)).*pGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform final calcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill scalling matrix
D = sparse(4*(numIntervals),4*(numIntervals));

for index = 3:4:length(D)-2
   D(index,index) = 1;
   D(index+1,index+1) = 10^4*numIntervals;
   D(index+2,index+2) = 1;
   D(index+3,index+3) = 10^4*numIntervals;
end
D(1,1)=1;
D(2,2)=1;
D(length(D)-1,length(D)-1)=1;
D(length(D),length(D))=1;

%Scale Matrices using scaling matrix
LHSscaled = D*massmatrix*D;
RHSscaled = D*massrhs;

%Solve simultaneous equations
ScaledCoeffs = LHSscaled\RHSscaled;
Coeffs = D'*ScaledCoeffs;

%Create second scaling matrix

%Collect Cj's in an array for display purposes
Cnjs = zeros(1,numIntervals+1);
Dnjs = zeros(1,numIntervals+1);
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);

for index=2:length(Cnjs)-1
   Cnjs(1,index) = Coeffs(4*index-5);
   Dnjs(1,index) = Coeffs(4*index-4);
   Cpjs(1,index) = Coeffs(4*index-3);
   Dpjs(1,index) = Coeffs(4*index-2);
end
Cnjs(1) = Coeffs(1);
Dnjs(1) = Coeffs(1)*alphan;
Cpjs(1) = Coeffs(2);
Dpjs(1) = Coeffs(2)*alphap;

Cnjs(numIntervals+1)=Coeffs(length(Coeffs)-1);
Dnjs(numIntervals+1)=Coeffs(length(Coeffs)-1)*betan;
Cpjs(numIntervals+1)=Coeffs(length(Coeffs));
Dpjs(numIntervals+1)=Coeffs(length(Coeffs))*betap;

InitialCnjs = Cnjs;
InitialDnjs = Dnjs;
InitialCpjs = Cpjs;
InitialDpjs = Dpjs;

nOrig=nGuessOrig;
pOrig=pGuessOrig;
trueCnjs = nOrig(node);
trueDnjs = nprimeOrig(node);
trueCpjs = pOrig(node);
trueDpjs = pprimeOrig(node);

fig4 = figure;
plot(indexn,Cnjs,'x',node,nGuess(node),indexn,nOrig(indexn)),title('nGuess (recomposition)'),legend('Cnjs','Error','Orginial');
saveas(fig4,strcat(filepath,'/nGuess_recomposition.jpg'))

fig4 = figure;
plot(indexn,Dnjs,'x',node,nprime(node),indexn,nprimeOrig(indexn)),title('nprime (recomposition)'),legend('Dnjs','Error','Orginial');
saveas(fig4,strcat(filepath,'/nprime_recomposition.jpg'))

fig5 = figure;
plot(indexn,Cpjs,'x',node,pGuess(node),indexn,pOrig(indexn)),title('pGuess (recomposition)'),legend('Cpjs','Error','Orginial');
saveas(fig5,strcat(filepath,'/pGuess_recomposition.jpg'))

fig5 = figure;
plot(indexn,Dpjs,'x',node,pprime(node),indexn,pprimeOrig(indexn)),title('pprime (recomposition)'),legend('Dpjs','Error','Orginial');
saveas(fig5,strcat(filepath,'/pprime_recomposition.jpg'))

%[output6,Jacobian6,RHSOut,Rn1,Rn2,Rn3,Rp1,Rp2,Rp3]=callFunction5b(Coeffs);

options = optimoptions('fsolve','Display','iter','TolFun',10^10,'Jacobian','on','TypicalX',Coeffs);
[x,fval,exitflag,output,Jacobian]= fsolve(@(x)callFunction5(x),Coeffs,options);

%Displays runtime of m-file
time = cputime-t



%Store variable values to be used later
save(strcat(filepath,'/',num2str(numIntervals),'terms',num2str(maxIterations),'iterations.mat'))