function [output,Jacobian] = callFunction4(inputCoeffs)

%Copy various constants and parameters used by body of code
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

startPoint = 0;
endPoint = Wn;
totalLength = endPoint-startPoint;

%Set number of intervals to divide domain into
numIntervals = 100;
intervalLength = totalLength/numIntervals;
dx = intervalLength;
% Store node locations
node = startPoint:dx:totalLength;

AbsError=10^-18;
RelError=10^-6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Gopp, Gopn
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


Gopn = @(x) -(cn1.*(n1*(gamman1*gamman1)*exp(-gamman1*x) + n2*(gamman2*gamman2)*exp(-gamman2*x) + n3*(gamman3*gamman3)*exp(-gamman3*x)+ n4*(gamman4*gamman4)*exp(-gamman4*x))    + cn2.*pprimeOrig(x)./pGuessOrig(x).*(-n1*(gamman1)*exp(-gamman1*x) - n2*(gamman2)*exp(-gamman2*x) - n3*(gamman3)*exp(-gamman3*x)- n4*(gamman4)*exp(-gamman4*x)));%+(ksiN*T*mun/q*(pprimeprimeOrig(x)./pGuessOrig(x) - (pprimeOrig(x)./pGuessOrig(x)).^2) - 1/taun).*(n1*exp(-gamman1*x) + n2*exp(-gamman2*x) + n3*exp(-gamman3*x));     
Gopp = @(x) -(cp1.*(p1*(gammap1*gammap1)*exp(-gammap1*x) + p2*(gammap2*gammap2)*exp(-gammap2*x) + p3*(gammap3*gammap3)*exp(-gammap3*x)+ p4*(gammap4*gammap4)*exp(-gammap4*x))    + cp2.*nprimeOrig(x)./nGuessOrig(x).*(-p1*(gammap1)*exp(-gammap1*x) - p2*(gammap2)*exp(-gammap2*x) - p3*(gammap3)*exp(-gammap3*x)- p4*(gammap4)*exp(-gammap4*x)));%+(ksiP*T*mup/q*((nprimeOrig(x)./nGuessOrig(x)).^2 - nprimeprimeOrig(x)./nGuessOrig(x)) - 1/taup).*(p1*exp(-gammap1*x) + p2*exp(-gammap2*x) + p3*exp(-gammap3*x));     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End Calculation of Gopp,Gopn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define general basis functions for interval 0 to 1
% Use inline functions; define generic basis functions and derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End Basis Function Defintions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Front end coefficients
Cn1=Dn - ksiN*T*mun/q;
Cn2=ksiN*T*mun/q;
Cn3=ksiN*T*mun/q;
Cp1=ksiP*T*mup/q - Dp;
Cp2=-ksiP*T*mup/q;
Cp3=ksiP*T*mup/q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphan=nprimeOrig(0)/nGuessOrig(0);
alphap=pprimeOrig(0)/pGuessOrig(0);
betan=nprimeOrig(Wn)/nGuessOrig(Wn);
betap=pprimeOrig(Wn)/pGuessOrig(Wn);

%Collect Cj's in an array for display purposes
Cnjs = zeros(1,numIntervals+1);
Dnjs = zeros(1,numIntervals+1);
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);

for index=2:length(Cnjs)-1
   Cnjs(1,index) = inputCoeffs(4*index-5);
   Dnjs(1,index) = inputCoeffs(4*index-4);
   Cpjs(1,index) = inputCoeffs(4*index-3);
   Dpjs(1,index) = inputCoeffs(4*index-2);
end

Cnjs(1) = inputCoeffs(1);
Dnjs(1) = inputCoeffs(1)*alphan;
Cpjs(1) = inputCoeffs(2);
Dpjs(1) = inputCoeffs(2)*alphap;

Cnjs(numIntervals+1)=inputCoeffs(length(inputCoeffs)-1);
Dnjs(numIntervals+1)=inputCoeffs(length(inputCoeffs)-1)*betan;
Cpjs(numIntervals+1)=inputCoeffs(length(inputCoeffs));
Dpjs(numIntervals+1)=inputCoeffs(length(inputCoeffs))*betap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates for RHS of equation (forcing function times test functions)
%This calculation is performed outside the iterative loop, as it will not
%change over the course of changing p(x), n(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs_matrix = zeros(4*(numIntervals),1);

for k = 1:numIntervals-2
    rhs_matrix(4*k-1) =      rhs_matrix(4*k-1) +        integral(@(x) f1((x-node(k+1))).*Gopn(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	rhs_matrix(4*k) =        rhs_matrix(4*k)   +        integral(@(x) f3((x-node(k+1))).*Gopn(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    rhs_matrix(4*k+3) =      rhs_matrix(4*k+3) +        integral(@(x) f2((x-node(k+1))).*Gopn(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	rhs_matrix(4*k+4) =      rhs_matrix(4*k+4) +        integral(@(x) f4((x-node(k+1))).*Gopn(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %pGuess terms
	rhs_matrix(4*k+1) =    rhs_matrix(4*k+1) +      integral(@(x) f1((x-node(k+1))).*Gopp(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	rhs_matrix(4*k+2) =    rhs_matrix(4*k+2) +      integral(@(x) f3((x-node(k+1))).*Gopp(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    rhs_matrix(4*k+5) =    rhs_matrix(4*k+5) +      integral(@(x) f2((x-node(k+1))).*Gopp(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	rhs_matrix(4*k+6) =    rhs_matrix(4*k+6) +      integral(@(x) f4((x-node(k+1))).*Gopp(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
       
end

%RHS nGuess edge cases
rhs_matrix(1) = rhs_matrix(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*Gopn(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(3) = rhs_matrix(3) + integral(@(x) f2((x-node(1))).*Gopn(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(4) = rhs_matrix(4) + integral(@(x) f4((x-node(1))).*Gopn(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%RHS nGuess edge cases
rhs_matrix(2) = rhs_matrix(2) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*Gopp(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(5) = rhs_matrix(5) + integral(@(x) f2((x-node(1))).*Gopp(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(6) = rhs_matrix(6) + integral(@(x) f4((x-node(1))).*Gopp(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

%RHS edge cases
%nGuess
rhs_matrix(last-1) = rhs_matrix(last-1)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*Gopn(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last-5) = rhs_matrix(last-5) + integral(@(x) f1((x-node(end)+dx)).*Gopn(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last-4) = rhs_matrix(last-4) + integral(@(x) f3((x-node(end)+dx)).*Gopn(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%pGuess
rhs_matrix(last) = rhs_matrix(last)+ integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*Gopp(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last-3) = rhs_matrix(last-3) + integral(@(x) f1((x-node(end)+dx)).*Gopp(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last-2) = rhs_matrix(last-2) + integral(@(x) f3((x-node(end)+dx)).*Gopp(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

RHSOut =rhs_matrix;

matrixDelNpp = sparse(4*(numIntervals),4*(numIntervals));

for k = 1:numIntervals-2
   k=k+1;
    %Cn1
    nBlock = @(x) Cn1;
    %-Cp3*pk/nk  
    pBlock = @(x) 0;%-Cp3.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
    k=k-1;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 1  - (Upper left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNpp(phi(n,r),xx) - nEq
	matrixDelNpp(4*k-1,4*k-1) = matrixDelNpp(4*k-1,4*k-1) + integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    matrixDelNpp(4*k,4*k-1)   = matrixDelNpp(4*k,4*k-1) +   integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(n,r),xx) - nEq
    matrixDelNpp(4*k-1,4*k) = matrixDelNpp(4*k-1,4*k) + integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k,4*k)   = matrixDelNpp(4*k,4*k) +   integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    

    %matrixDelNpp(phi(p,r),xx) - pEq
	matrixDelNpp(4*k+1,4*k-1) = matrixDelNpp(4*k+1,4*k-1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+2,4*k-1) = matrixDelNpp(4*k+2,4*k-1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(p,r),xx) - pEq
	matrixDelNpp(4*k+1,4*k) = matrixDelNpp(4*k+1,4*k) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+2,4*k) = matrixDelNpp(4*k+2,4*k) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2  - (Upper right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNpp(phi(n,r),xx) - nEq
	matrixDelNpp(4*k-1,4*k+3) = matrixDelNpp(4*k-1,4*k+3) + integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k,4*k+3)   = matrixDelNpp(4*k,4*k+3) +   integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(n,r),xx) - nEq
    matrixDelNpp(4*k-1,4*k+4) = matrixDelNpp(4*k-1,4*k+4) + integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k,4*k+4)   = matrixDelNpp(4*k,4*k+4) +   integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(phi(p,r),xx) - pEq
	matrixDelNpp(4*k+1,4*k+3) = matrixDelNpp(4*k+1,4*k+3) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+2,4*k+3) = matrixDelNpp(4*k+2,4*k+3) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(p,r),xx) - pEq
	matrixDelNpp(4*k+1,4*k+4) = matrixDelNpp(4*k+1,4*k+4) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+2,4*k+4) = matrixDelNpp(4*k+2,4*k+4) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 3  - (Bottom left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNpp(phi(n,r),xx) - nEq
	matrixDelNpp(4*k+3,4*k-1) = matrixDelNpp(4*k+3,4*k-1) + integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+4,4*k-1) = matrixDelNpp(4*k+4,4*k-1) + integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(n,r),xx) - nEq
    matrixDelNpp(4*k+3,4*k) = matrixDelNpp(4*k+3,4*k) + integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+4,4*k) = matrixDelNpp(4*k+4,4*k) + integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(phi(p,r),xx) - pEq
	matrixDelNpp(4*k+5,4*k-1) = matrixDelNpp(4*k+5,4*k-1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+6,4*k-1) = matrixDelNpp(4*k+6,4*k-1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(p,r),xx) - pEq
	matrixDelNpp(4*k+5,4*k) = matrixDelNpp(4*k+5,4*k) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+6,4*k) = matrixDelNpp(4*k+6,4*k) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 4  - (Bottom right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNpp(phi(n,r),xx) - nEq
	matrixDelNpp(4*k+3,4*k+3) = matrixDelNpp(4*k+3,4*k+3) + integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    matrixDelNpp(4*k+4,4*k+3) = matrixDelNpp(4*k+4,4*k+3) + integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(n,r),xx) - nEq
    matrixDelNpp(4*k+3,4*k+4) = matrixDelNpp(4*k+3,4*k+4) + integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+4,4*k+4) = matrixDelNpp(4*k+4,4*k+4) + integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(phi(p,r),xx) - pEq
	matrixDelNpp(4*k+5,4*k+3) = matrixDelNpp(4*k+5,4*k+3) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+6,4*k+3) = matrixDelNpp(4*k+6,4*k+3) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNpp(psi(p,r),xx) - pEq
	matrixDelNpp(4*k+5,4*k+4) = matrixDelNpp(4*k+5,4*k+4) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNpp(4*k+6,4*k+4) = matrixDelNpp(4*k+6,4*k+4) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
end    
    k=1;
    nBlock = @(x) Cn1;
    pBlock = @(x) 0;%-Cp3.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
      
%Handle 1st interval 'half case'
%Top left 2x2 block
matrixDelNpp(1,1)= integral(@(x) nBlock(x).*(f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(2,1)= integral(@(x) pBlock(x).*(f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Left 2x4 block
matrixDelNpp(3,1) = integral(@(x) nBlock(x).*(f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(4,1) = integral(@(x) nBlock(x).*(f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelNpp(5,1) = integral(@(x) pBlock(x).*(f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(6,1) = integral(@(x) pBlock(x).*(f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Right 4x2 block
matrixDelNpp(1,3) = integral(@(x) nBlock(x).*f2primeprime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(1,4) = integral(@(x) nBlock(x).*f4primeprime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelNpp(2,3) = integral(@(x) pBlock(x).*f2primeprime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(2,4) = integral(@(x) pBlock(x).*f4primeprime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelNpp(phi(n,r),xx)
matrixDelNpp(3,3) = matrixDelNpp(3,3) + integral(@(x) nBlock(x).*f2primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(4,3) = matrixDelNpp(4,3) + integral(@(x) nBlock(x).*f2primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNpp(psi(n,r),xx)
matrixDelNpp(3,4) = matrixDelNpp(3,4) + integral(@(x) nBlock(x).*f4primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(4,4) = matrixDelNpp(4,4) + integral(@(x) nBlock(x).*f4primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNpp(phi(p,r),xx)
matrixDelNpp(5,3) = matrixDelNpp(5,3) + integral(@(x) pBlock(x).*f2primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(6,3) = matrixDelNpp(6,3) + integral(@(x) pBlock(x).*f2primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNpp(psi(p,r),xx)
matrixDelNpp(5,4) = matrixDelNpp(5,4) + integral(@(x) pBlock(x).*f4primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(6,4) = matrixDelNpp(6,4) + integral(@(x) pBlock(x).*f4primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

    k=numIntervals;
    nBlock = @(x) Cn1;
    pBlock = @(x) 0;%-Cp3.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
    
%Handle 1st interval 'half case'
%Bottom right 2x2 block
matrixDelNpp(last,last-1) =   integral(@(x) pBlock(x).*(f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-1,last-1) = integral(@(x) nBlock(x).*(f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Bottom 2x4 block
matrixDelNpp(last,last-4) = integral(@(x) pBlock(x).*f3primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last,last-5) = integral(@(x) pBlock(x).*f1primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelNpp(last-1,last-4) = integral(@(x) nBlock(x).*f3primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-1,last-5) = integral(@(x) nBlock(x).*f1primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Left 4x2 Block
matrixDelNpp(last-2,last-1) = integral(@(x) pBlock(x).*(f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-3,last-1) = integral(@(x) pBlock(x).*(f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelNpp(last-4,last-1) = integral(@(x) nBlock(x).*(f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-5,last-1) = integral(@(x) nBlock(x).*(f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelNpp(phi(n,r),xx)
matrixDelNpp(last-5,last-5) = matrixDelNpp(last-5,last-5) + integral(@(x) nBlock(x).*f1primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-4,last-5) = matrixDelNpp(last-4,last-5) + integral(@(x) nBlock(x).*f1primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNpp(psi(n,r),xx)
matrixDelNpp(last-5,last-4) = matrixDelNpp(last-5,last-4) + integral(@(x) nBlock(x).*f3primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-4,last-4) = matrixDelNpp(last-4,last-4) + integral(@(x) nBlock(x).*f3primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNpp(phi(p,r),xx)
matrixDelNpp(last-3,last-5) = matrixDelNpp(last-3,last-5) + integral(@(x) pBlock(x).*f1primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-2,last-5) = matrixDelNpp(last-2,last-5) + integral(@(x) pBlock(x).*f1primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNpp(psi(p,r),xx)
matrixDelNpp(last-3,last-4) = matrixDelNpp(last-3,last-4) + integral(@(x) pBlock(x).*f3primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNpp(last-2,last-4) = matrixDelNpp(last-2,last-4) + integral(@(x) pBlock(x).*f3primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
   
Matrix1Finaln = matrixDelNpp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End RHS calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with (delp)''(x) term
% %Performed outside iterative loop because this term is not influenced by
% %p(x) or n(x).  Should not change over calculations.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Front end coefficient of n''(x) term

matrixDelPpp = sparse(4*(numIntervals),4*(numIntervals));

for k = 1:numIntervals-2
    k=k+1;
    %Cn3*nk/pk  
    nBlock = @(x) 0;%Cn3.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    %Cp1 
    pBlock = @(x) Cp1;
    k=k-1;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 1  - (Upper left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPpp(phi(n,r),xx) - nEq
	matrixDelPpp(4*k-1,4*k+1) = matrixDelPpp(4*k-1,4*k+1) + integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k,4*k+1)   = matrixDelPpp(4*k,4*k+1) +   integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	
    %matrixDelPpp(psi(n,r),xx) - nEq
    matrixDelPpp(4*k-1,4*k+2) = matrixDelPpp(4*k-1,4*k+2) + integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k,4*k+2)   = matrixDelPpp(4*k,4*k+2) +   integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	
    %matrixDelPpp(phi(p,r),xx) - pEq
    matrixDelPpp(4*k+1,4*k+1) = matrixDelPpp(4*k+1,4*k+1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    matrixDelPpp(4*k+2,4*k+1) = matrixDelPpp(4*k+2,4*k+1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(psi(p,r),xx) - pEq
    matrixDelPpp(4*k+1,4*k+2) = matrixDelPpp(4*k+1,4*k+2) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+2,4*k+2) = matrixDelPpp(4*k+2,4*k+2) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2  - (Upper right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPpp(phi(n,r),xx) - nEq
	matrixDelPpp(4*k-1,4*k+5) = matrixDelPpp(4*k-1,4*k+5) + integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k,4*k+5)   = matrixDelPpp(4*k,4*k+5) +   integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	
    %matrixDelPpp(psi(n,r),xx) - nEq
    matrixDelPpp(4*k-1,4*k+6) = matrixDelPpp(4*k-1,4*k+6) + integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k,4*k+6)   = matrixDelPpp(4*k,4*k+6) +   integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	
    %matrixDelPpp(phi(p,r),xx) - pEq
    matrixDelPpp(4*k+1,4*k+5) = matrixDelPpp(4*k+1,4*k+5) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+2,4*k+5) = matrixDelPpp(4*k+2,4*k+5) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(psi(p,r),xx) - pEq
    matrixDelPpp(4*k+1,4*k+6) = matrixDelPpp(4*k+1,4*k+6) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+2,4*k+6) = matrixDelPpp(4*k+2,4*k+6) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 3  - (Bottom left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPpp(phi(n,r),xx) - nEq
	matrixDelPpp(4*k+3,4*k+1) = matrixDelPpp(4*k+3,4*k+1) + integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+4,4*k+1) = matrixDelPpp(4*k+4,4*k+1) + integral(@(x) nBlock(x).*f1primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(psi(n,r),xx) - nEq
    matrixDelPpp(4*k+3,4*k+2) = matrixDelPpp(4*k+3,4*k+2) + integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+4,4*k+2) = matrixDelPpp(4*k+4,4*k+2) + integral(@(x) nBlock(x).*f3primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(phi(p,r),xx) - pEq
	matrixDelPpp(4*k+5,4*k+1) = matrixDelPpp(4*k+5,4*k+1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+6,4*k+1) = matrixDelPpp(4*k+6,4*k+1) + integral(@(x) pBlock(x).*f1primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(psi(p,r),xx) - pEq
	matrixDelPpp(4*k+5,4*k+2) = matrixDelPpp(4*k+5,4*k+2) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+6,4*k+2) = matrixDelPpp(4*k+6,4*k+2) + integral(@(x) pBlock(x).*f3primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 4  - (Bottom right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPpp(phi(n,r),xx) - nEq
	matrixDelPpp(4*k+3,4*k+5) = matrixDelPpp(4*k+3,4*k+5) + integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+4,4*k+5) = matrixDelPpp(4*k+4,4*k+5) + integral(@(x) nBlock(x).*f2primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(psi(n,r),xx) - nEq
    matrixDelPpp(4*k+3,4*k+6) = matrixDelPpp(4*k+3,4*k+6) + integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+4,4*k+6) = matrixDelPpp(4*k+4,4*k+6) + integral(@(x) nBlock(x).*f4primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(phi(p,r),xx) - pEq
	matrixDelPpp(4*k+5,4*k+5) = matrixDelPpp(4*k+5,4*k+5) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    matrixDelPpp(4*k+6,4*k+5) = matrixDelPpp(4*k+6,4*k+5) + integral(@(x) pBlock(x).*f2primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPpp(psi(p,r),xx) - pEq
    matrixDelPpp(4*k+5,4*k+6) = matrixDelPpp(4*k+5,4*k+6) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPpp(4*k+6,4*k+6) = matrixDelPpp(4*k+6,4*k+6) + integral(@(x) pBlock(x).*f4primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
end    
    k=1;
    nBlock = @(x) 0;%Cn3.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) Cp1;
    
%Handle 1st interval 'half case'
%Top left 2x2 block
matrixDelPpp(1,2)= integral(@(x) nBlock(x).*(f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(2,2)= integral(@(x) pBlock(x).*(f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Left 2x4 block
matrixDelPpp(3,2) = integral(@(x) nBlock(x).*(f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(4,2) = integral(@(x) nBlock(x).*(f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelPpp(5,2) = integral(@(x) pBlock(x).*(f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(6,2) = integral(@(x) pBlock(x).*(f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Right 4x2 block
matrixDelPpp(1,5) = integral(@(x) nBlock(x).*f2primeprime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(1,6) = integral(@(x) nBlock(x).*f4primeprime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelPpp(2,5) = integral(@(x) pBlock(x).*f2primeprime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(2,6) = integral(@(x) pBlock(x).*f4primeprime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelPpp(phi(n,r),xx)
matrixDelPpp(3,5) = matrixDelPpp(3,5) + integral(@(x) nBlock(x).*f2primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(4,5) = matrixDelPpp(4,5) + integral(@(x) nBlock(x).*f2primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPpp(psi(n,r),xx)
matrixDelPpp(3,6) = matrixDelPpp(3,6) + integral(@(x) nBlock(x).*f4primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(4,6) = matrixDelPpp(4,6) + integral(@(x) nBlock(x).*f4primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPpp(phi(p,r),xx)
matrixDelPpp(5,5) = matrixDelPpp(5,5) + integral(@(x) pBlock(x).*f2primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(6,5) = matrixDelPpp(6,5) + integral(@(x) pBlock(x).*f2primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPpp(psi(p,r),xx)
matrixDelPpp(5,6) = matrixDelPpp(5,6) + integral(@(x) pBlock(x).*f4primeprime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(6,6) = matrixDelPpp(6,6) + integral(@(x) pBlock(x).*f4primeprime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

    k=numIntervals;
    nBlock = @(x) 0;%Cn3.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) Cp1;
    
%Handle 1st interval 'half case'
%Bottom right 2x2 block
matrixDelPpp(last,last) =   integral(@(x) pBlock(x).*(f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-1,last) = integral(@(x) nBlock(x).*(f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Bottom 2x4 block
matrixDelPpp(last,last-2) = integral(@(x) pBlock(x).*f3primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last,last-3) = integral(@(x) pBlock(x).*f1primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelPpp(last-1,last-2) = integral(@(x) nBlock(x).*f3primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-1,last-3) = integral(@(x) nBlock(x).*f1primeprime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Right 4x2 Block
matrixDelPpp(last-2,last) = integral(@(x) pBlock(x).*(f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-3,last) = integral(@(x) pBlock(x).*(f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelPpp(last-4,last) = integral(@(x) nBlock(x).*(f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-5,last) = integral(@(x) nBlock(x).*(f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelPpp(phi(n,r),xx)
matrixDelPpp(last-5,last-3) = matrixDelPpp(last-5,last-3) + integral(@(x) nBlock(x).*f1primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-4,last-3) = matrixDelPpp(last-4,last-3) + integral(@(x) nBlock(x).*f1primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPpp(psi(n,r),xx)
matrixDelPpp(last-5,last-2) = matrixDelPpp(last-5,last-2) + integral(@(x) nBlock(x).*f3primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-4,last-2) = matrixDelPpp(last-4,last-2) + integral(@(x) nBlock(x).*f3primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPpp(phi(p,r),xx)
matrixDelPpp(last-3,last-3) = matrixDelPpp(last-3,last-3) + integral(@(x) pBlock(x).*f1primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-2,last-3) = matrixDelPpp(last-2,last-3) + integral(@(x) pBlock(x).*f1primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPpp(psi(p,r),xx)
matrixDelPpp(last-3,last-2) = matrixDelPpp(last-3,last-2) + integral(@(x) pBlock(x).*f3primeprime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPpp(last-2,last-2) = matrixDelPpp(last-2,last-2) + integral(@(x) pBlock(x).*f3primeprime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
   
Matrix1Finalp = matrixDelPpp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of (delp)''(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with (deln)'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Front end coefficient for this matrix.  DOES NOT INCLUDE CONTRIBUTION FROM
%p(x) DISTRIBUTION.

   matrixDelNp = sparse(4*(numIntervals),4*(numIntervals));


for k = 1:numIntervals-2
   k=k+1;
    %Cn2*p'k/pk
    nBlock = @(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    %Cp2 p'k/nk + 2*Cn3*(n'k) pk/(nk)^2   
    pBlock = @(x) Cp2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
    k=k-1;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 1  - (Upper left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNp(phi(n,r),xx) - nEq
	matrixDelNp(4*k-1,4*k-1) = matrixDelNp(4*k-1,4*k-1) + integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k,4*k-1)   = matrixDelNp(4*k,4*k-1) +   integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(n,r),xx) - nEq
    matrixDelNp(4*k-1,4*k) = matrixDelNp(4*k-1,4*k) + integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k,4*k)   = matrixDelNp(4*k,4*k) +   integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(phi(p,r),xx) - pEq
	matrixDelNp(4*k+1,4*k-1) = matrixDelNp(4*k+1,4*k-1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+2,4*k-1) = matrixDelNp(4*k+2,4*k-1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(p,r),xx) - pEq
	matrixDelNp(4*k+1,4*k) = matrixDelNp(4*k+1,4*k) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+2,4*k) = matrixDelNp(4*k+2,4*k) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2  - (Upper right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNp(phi(n,r),xx) - nEq
	matrixDelNp(4*k-1,4*k+3) = matrixDelNp(4*k-1,4*k+3) + integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k,4*k+3)   = matrixDelNp(4*k,4*k+3) +   integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(n,r),xx) - nEq
    matrixDelNp(4*k-1,4*k+4) = matrixDelNp(4*k-1,4*k+4) + integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k,4*k+4)   = matrixDelNp(4*k,4*k+4) +   integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(phi(p,r),xx) - pEq
	matrixDelNp(4*k+1,4*k+3) = matrixDelNp(4*k+1,4*k+3) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+2,4*k+3) = matrixDelNp(4*k+2,4*k+3) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(p,r),xx) - pEq
	matrixDelNp(4*k+1,4*k+4) = matrixDelNp(4*k+1,4*k+4) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+2,4*k+4) = matrixDelNp(4*k+2,4*k+4) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 3  - (Bottom left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNp(phi(n,r),xx) - nEq
	matrixDelNp(4*k+3,4*k-1) = matrixDelNp(4*k+3,4*k-1) + integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+4,4*k-1) = matrixDelNp(4*k+4,4*k-1) + integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(n,r),xx) - nEq
    matrixDelNp(4*k+3,4*k) = matrixDelNp(4*k+3,4*k) + integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+4,4*k) = matrixDelNp(4*k+4,4*k) + integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(phi(p,r),xx) - pEq
	matrixDelNp(4*k+5,4*k-1) = matrixDelNp(4*k+5,4*k-1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+6,4*k-1) = matrixDelNp(4*k+6,4*k-1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(p,r),xx) - pEq
	matrixDelNp(4*k+5,4*k) = matrixDelNp(4*k+5,4*k) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+6,4*k) = matrixDelNp(4*k+6,4*k) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 4  - (Bottom right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelNp(phi(n,r),xx) - nEq
	matrixDelNp(4*k+3,4*k+3) = matrixDelNp(4*k+3,4*k+3) + integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+4,4*k+3) = matrixDelNp(4*k+4,4*k+3) + integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(n,r),xx) - nEq
    matrixDelNp(4*k+3,4*k+4) = matrixDelNp(4*k+3,4*k+4) + integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+4,4*k+4) = matrixDelNp(4*k+4,4*k+4) + integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(phi(p,r),xx) - pEq
	matrixDelNp(4*k+5,4*k+3) = matrixDelNp(4*k+5,4*k+3) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+6,4*k+3) = matrixDelNp(4*k+6,4*k+3) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelNp(psi(p,r),xx) - pEq
	matrixDelNp(4*k+5,4*k+4) = matrixDelNp(4*k+5,4*k+4) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelNp(4*k+6,4*k+4) = matrixDelNp(4*k+6,4*k+4) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
end    
    k=1;
    nBlock = @(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) Cp2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
        
%Handle 1st interval 'half case'
%Top left 2x2 block
matrixDelNp(1,1)= integral(@(x) nBlock(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(2,1)= integral(@(x) pBlock(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Left 2x4 block
matrixDelNp(3,1) = integral(@(x) nBlock(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(4,1) = integral(@(x) nBlock(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelNp(5,1) = integral(@(x) pBlock(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(6,1) = integral(@(x) pBlock(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Right 4x2 block
matrixDelNp(1,3) = integral(@(x) nBlock(x).*f2prime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(1,4) = integral(@(x) nBlock(x).*f4prime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelNp(2,3) = integral(@(x) pBlock(x).*f2prime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(2,4) = integral(@(x) pBlock(x).*f4prime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelNp(phi(n,r),xx)
matrixDelNp(3,3) = matrixDelNp(3,3) + integral(@(x) nBlock(x).*f2prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(4,3) = matrixDelNp(4,3) + integral(@(x) nBlock(x).*f2prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNp(psi(n,r),xx)
matrixDelNp(3,4) = matrixDelNp(3,4) + integral(@(x) nBlock(x).*f4prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(4,4) = matrixDelNp(4,4) + integral(@(x) nBlock(x).*f4prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNp(phi(p,r),xx)
matrixDelNp(5,3) = matrixDelNp(5,3) + integral(@(x) pBlock(x).*f2prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(6,3) = matrixDelNp(6,3) + integral(@(x) pBlock(x).*f2prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNp(psi(p,r),xx)
matrixDelNp(5,4) = matrixDelNp(5,4) + integral(@(x) pBlock(x).*f4prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(6,4) = matrixDelNp(6,4) + integral(@(x) pBlock(x).*f4prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

    k=numIntervals;
    nBlock = @(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) Cp2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
      
%Handle 1st interval 'half case'
%Bottom right 2x2 block
matrixDelNp(last,last-1) =   integral(@(x) pBlock(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-1,last-1) = integral(@(x) nBlock(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Bottom 2x4 block
matrixDelNp(last,last-4) = integral(@(x) pBlock(x).*f3prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last,last-5) = integral(@(x) pBlock(x).*f1prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelNp(last-1,last-4) = integral(@(x) nBlock(x).*f3prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-1,last-5) = integral(@(x) nBlock(x).*f1prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Left 4x2 Block
matrixDelNp(last-2,last-1) = integral(@(x) pBlock(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-3,last-1) = integral(@(x) pBlock(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelNp(last-4,last-1) = integral(@(x) nBlock(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-5,last-1) = integral(@(x) nBlock(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelNp(phi(n,r),xx)
matrixDelNp(last-5,last-5) = matrixDelNp(last-5,last-5) + integral(@(x) nBlock(x).*f1prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-4,last-5) = matrixDelNp(last-4,last-5) + integral(@(x) nBlock(x).*f1prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNp(psi(n,r),xx)
matrixDelNp(last-5,last-4) = matrixDelNp(last-5,last-4) + integral(@(x) nBlock(x).*f3prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-4,last-4) = matrixDelNp(last-4,last-4) + integral(@(x) nBlock(x).*f3prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNp(phi(p,r),xx)
matrixDelNp(last-3,last-5) = matrixDelNp(last-3,last-5) + integral(@(x) pBlock(x).*f1prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-2,last-5) = matrixDelNp(last-2,last-5) + integral(@(x) pBlock(x).*f1prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelNp(psi(p,r),xx)
matrixDelNp(last-3,last-4) = matrixDelNp(last-3,last-4) + integral(@(x) pBlock(x).*f3prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelNp(last-2,last-4) = matrixDelNp(last-2,last-4) + integral(@(x) pBlock(x).*f3prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
   

   Matrix2Finaln = matrixDelNp;
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of deln'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with delp'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Front end coefficient for this matrix.

matrixDelPp = sparse(4*(numIntervals),4*(numIntervals));

for k = 1:numIntervals-2
   k=k+1;
    %Cn2*nk'/pk - 2*Cn3*pk'*nk/(pk)^2
    nBlock = @(x) Cn2.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))/(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    %Cp2*nk'/nk   
    pBlock = @(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
    k=k-1;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 1  - (Upper left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPp(phi(n,r),xx) - nEq
	matrixDelPp(4*k-1,4*k+1) = matrixDelPp(4*k-1,4*k+1) + integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k,4*k+1)   = matrixDelPp(4*k,4*k+1) +   integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(n,r),xx) - nEq
    matrixDelPp(4*k-1,4*k+2) = matrixDelPp(4*k-1,4*k+2) + integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k,4*k+2)   = matrixDelPp(4*k,4*k+2) +   integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(phi(p,r),xx) - pEq
	matrixDelPp(4*k+1,4*k+1) = matrixDelPp(4*k+1,4*k+1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+2,4*k+1) = matrixDelPp(4*k+2,4*k+1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(p,r),xx) - pEq
	matrixDelPp(4*k+1,4*k+2) = matrixDelPp(4*k+1,4*k+2) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+2,4*k+2) = matrixDelPp(4*k+2,4*k+2) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2  - (Upper right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPp(phi(n,r),xx) - nEq
	matrixDelPp(4*k-1,4*k+5) = matrixDelPp(4*k-1,4*k+5) + integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k,4*k+5)   = matrixDelPp(4*k,4*k+5) +   integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(n,r),xx) - nEq
    matrixDelPp(4*k-1,4*k+6) = matrixDelPp(4*k-1,4*k+6) + integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k,4*k+6)   = matrixDelPp(4*k,4*k+6) +   integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(phi(p,r),xx) - pEq
	matrixDelPp(4*k+1,4*k+5) = matrixDelPp(4*k+1,4*k+5) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+2,4*k+5) = matrixDelPp(4*k+2,4*k+5) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(p,r),xx) - pEq
	matrixDelPp(4*k+1,4*k+6) = matrixDelPp(4*k+1,4*k+6) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+2,4*k+6) = matrixDelPp(4*k+2,4*k+6) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 3  - (Bottom left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPp(phi(n,r),xx) - nEq
	matrixDelPp(4*k+3,4*k+1) = matrixDelPp(4*k+3,4*k+1) + integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+4,4*k+1) = matrixDelPp(4*k+4,4*k+1) + integral(@(x) nBlock(x).*f1prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(n,r),xx) - nEq
    matrixDelPp(4*k+3,4*k+2) = matrixDelPp(4*k+3,4*k+2) + integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+4,4*k+2) = matrixDelPp(4*k+4,4*k+2) + integral(@(x) nBlock(x).*f3prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(phi(p,r),xx) - pEq
	matrixDelPp(4*k+5,4*k+1) = matrixDelPp(4*k+5,4*k+1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+6,4*k+1) = matrixDelPp(4*k+6,4*k+1) + integral(@(x) pBlock(x).*f1prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(p,r),xx) - pEq
	matrixDelPp(4*k+5,4*k+2) = matrixDelPp(4*k+5,4*k+2) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+6,4*k+2) = matrixDelPp(4*k+6,4*k+2) + integral(@(x) pBlock(x).*f3prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 4  - (Bottom right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelPp(phi(n,r),xx) - nEq
	matrixDelPp(4*k+3,4*k+5) = matrixDelPp(4*k+3,4*k+5) + integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+4,4*k+5) = matrixDelPp(4*k+4,4*k+5) + integral(@(x) nBlock(x).*f2prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(n,r),xx) - nEq
    matrixDelPp(4*k+3,4*k+6) = matrixDelPp(4*k+3,4*k+6) + integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+4,4*k+6) = matrixDelPp(4*k+4,4*k+6) + integral(@(x) nBlock(x).*f4prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(phi(p,r),xx) - pEq
	matrixDelPp(4*k+5,4*k+5) = matrixDelPp(4*k+5,4*k+5) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+6,4*k+5) = matrixDelPp(4*k+6,4*k+5) + integral(@(x) pBlock(x).*f2prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelPp(psi(p,r),xx) - pEq
	matrixDelPp(4*k+5,4*k+6) = matrixDelPp(4*k+5,4*k+6) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelPp(4*k+6,4*k+6) = matrixDelPp(4*k+6,4*k+6) + integral(@(x) pBlock(x).*f4prime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
end    
    k=1;
    nBlock = @(x) Cn2.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))/(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
         
%Handle 1st interval 'half case'
%Top left 2x2 block
matrixDelPp(1,2)= integral(@(x) nBlock(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(2,2)= integral(@(x) pBlock(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Left 2x4 block
matrixDelPp(3,2) = integral(@(x) nBlock(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(4,2) = integral(@(x) nBlock(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelPp(5,2) = integral(@(x) pBlock(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(6,2) = integral(@(x) pBlock(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Right 4x2 block
matrixDelPp(1,5) = integral(@(x) nBlock(x).*f2prime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(1,6) = integral(@(x) nBlock(x).*f4prime((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelPp(2,5) = integral(@(x) pBlock(x).*f2prime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(2,6) = integral(@(x) pBlock(x).*f4prime((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelPp(phi(n,r),xx)
matrixDelPp(3,5) = matrixDelPp(3,5) + integral(@(x) nBlock(x).*f2prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(4,5) = matrixDelPp(4,5) + integral(@(x) nBlock(x).*f2prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPp(psi(n,r),xx)
matrixDelPp(3,6) = matrixDelPp(3,6) + integral(@(x) nBlock(x).*f4prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(4,6) = matrixDelPp(4,6) + integral(@(x) nBlock(x).*f4prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPp(phi(p,r),xx)
matrixDelPp(5,5) = matrixDelPp(5,5) + integral(@(x) pBlock(x).*f2prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(6,5) = matrixDelPp(6,5) + integral(@(x) pBlock(x).*f2prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPp(psi(p,r),xx)
matrixDelPp(5,6) = matrixDelPp(5,6) + integral(@(x) pBlock(x).*f4prime((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(6,6) = matrixDelPp(6,6) + integral(@(x) pBlock(x).*f4prime((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

    k=numIntervals;
    nBlock = @(x) Cn2.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))/(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
         
%Handle 1st interval 'half case'
%Bottom right 2x2 block
matrixDelPp(last,last) =   integral(@(x) pBlock(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-1,last) = integral(@(x) nBlock(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Bottom 2x4 block
matrixDelPp(last,last-2) = integral(@(x) pBlock(x).*f3prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last,last-3) = integral(@(x) pBlock(x).*f1prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelPp(last-1,last-2) = integral(@(x) nBlock(x).*f3prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-1,last-3) = integral(@(x) nBlock(x).*f1prime((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Left 4x2 Block
matrixDelPp(last-2,last) = integral(@(x) pBlock(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-3,last) = integral(@(x) pBlock(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelPp(last-4,last) = integral(@(x) nBlock(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-5,last) = integral(@(x) nBlock(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelPp(phi(n,r),xx)
matrixDelPp(last-5,last-3) = matrixDelPp(last-5,last-3) + integral(@(x) nBlock(x).*f1prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-4,last-3) = matrixDelPp(last-4,last-3) + integral(@(x) nBlock(x).*f1prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPp(psi(n,r),xx)
matrixDelPp(last-5,last-2) = matrixDelPp(last-5,last-2) + integral(@(x) nBlock(x).*f3prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-4,last-2) = matrixDelPp(last-4,last-2) + integral(@(x) nBlock(x).*f3prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPp(phi(p,r),xx)
matrixDelPp(last-3,last-3) = matrixDelPp(last-3,last-3) + integral(@(x) pBlock(x).*f1prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-2,last-3) = matrixDelPp(last-2,last-3) + integral(@(x) pBlock(x).*f1prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelPp(psi(p,r),xx)
matrixDelPp(last-3,last-2) = matrixDelPp(last-3,last-2) + integral(@(x) pBlock(x).*f3prime(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelPp(last-2,last-2) = matrixDelPp(last-2,last-2) + integral(@(x) pBlock(x).*f3prime(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
   
Matrix2Finalp = matrixDelPp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of delp'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with deln(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrixDelN = sparse(4*(numIntervals),4*(numIntervals));

for k = 1:numIntervals-2
   k=k+1;
    %Cn3*(pk''/pk-(pk'/pk)^2)-1/taun  
    nBlock = @(x) 0;%Cn3.*((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))) - ((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2)-1/taun;
    %Cp3*nk''pk/nk^2-Cp2*nk' pk'/nk^2 - 2Cp3(nk')^2*pk/(nk)^3
    pBlock = @(x) -Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;
    k=k-1;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 1  - (Upper left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelN(phi(n,r),xx) - nEq
	matrixDelN(4*k-1,4*k-1) = matrixDelN(4*k-1,4*k-1) + integral(@(x) nBlock(x).*f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k,4*k-1)   = matrixDelN(4*k,4*k-1) +   integral(@(x) nBlock(x).*f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(psi(n,r),xx) - nEq
    matrixDelN(4*k-1,4*k) = matrixDelN(4*k-1,4*k) + integral(@(x) nBlock(x).*f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k,4*k)   = matrixDelN(4*k,4*k) +   integral(@(x) nBlock(x).*f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(phi(p,r),xx) - pEq
	matrixDelN(4*k+1,4*k-1) = matrixDelN(4*k+1,4*k-1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+2,4*k-1) = matrixDelN(4*k+2,4*k-1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(psi(p,r),xx) - pEq
	matrixDelN(4*k+1,4*k) = matrixDelN(4*k+1,4*k) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+2,4*k) = matrixDelN(4*k+2,4*k) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2  - (Upper right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelN(phi(n,r),xx) - nEq
	matrixDelN(4*k-1,4*k+3) = matrixDelN(4*k-1,4*k+3) + integral(@(x) nBlock(x).*f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k,4*k+3)   = matrixDelN(4*k,4*k+3) +   integral(@(x) nBlock(x).*f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(psi(n,r),xx) - nEq
    matrixDelN(4*k-1,4*k+4) = matrixDelN(4*k-1,4*k+4) + integral(@(x) nBlock(x).*f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k,4*k+4)   = matrixDelN(4*k,4*k+4) +   integral(@(x) nBlock(x).*f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(phi(p,r),xx) - pEq
	matrixDelN(4*k+1,4*k+3) = matrixDelN(4*k+1,4*k+3) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+2,4*k+3) = matrixDelN(4*k+2,4*k+3) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(psi(p,r),xx) - pEq
	matrixDelN(4*k+1,4*k+4) = matrixDelN(4*k+1,4*k+4) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+2,4*k+4) = matrixDelN(4*k+2,4*k+4) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 3  - (Bottom left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelN(phi(n,r),xx) - nEq
	matrixDelN(4*k+3,4*k-1) = matrixDelN(4*k+3,4*k-1) + integral(@(x) nBlock(x).*f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+4,4*k-1) = matrixDelN(4*k+4,4*k-1) + integral(@(x) nBlock(x).*f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    %matrixDelN(psi(n,r),xx) - nEq
    matrixDelN(4*k+3,4*k) = matrixDelN(4*k+3,4*k) + integral(@(x) nBlock(x).*f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+4,4*k) = matrixDelN(4*k+4,4*k) + integral(@(x) nBlock(x).*f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(phi(p,r),xx) - pEq
	matrixDelN(4*k+5,4*k-1) = matrixDelN(4*k+5,4*k-1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+6,4*k-1) = matrixDelN(4*k+6,4*k-1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(psi(p,r),xx) - pEq
	matrixDelN(4*k+5,4*k) = matrixDelN(4*k+5,4*k) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+6,4*k) = matrixDelN(4*k+6,4*k) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 4  - (Bottom right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelN(phi(n,r),xx) - nEq
	matrixDelN(4*k+3,4*k+3) = matrixDelN(4*k+3,4*k+3) + integral(@(x) nBlock(x).*f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+4,4*k+3) = matrixDelN(4*k+4,4*k+3) + integral(@(x) nBlock(x).*f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(psi(n,r),xx) - nEq
    matrixDelN(4*k+3,4*k+4) = matrixDelN(4*k+3,4*k+4) + integral(@(x) nBlock(x).*f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+4,4*k+4) = matrixDelN(4*k+4,4*k+4) + integral(@(x) nBlock(x).*f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(phi(p,r),xx) - pEq
	matrixDelN(4*k+5,4*k+3) = matrixDelN(4*k+5,4*k+3) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+6,4*k+3) = matrixDelN(4*k+6,4*k+3) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelN(psi(p,r),xx) - pEq
	matrixDelN(4*k+5,4*k+4) = matrixDelN(4*k+5,4*k+4) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelN(4*k+6,4*k+4) = matrixDelN(4*k+6,4*k+4) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
end    
    k=1;
    nBlock = @(x) 0;%Cn3.*((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))) - ((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2)-1/taun;
    pBlock = @(x) -Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;
    
%Handle 1st interval 'half case'
%Top left 2x2 block
matrixDelN(1,1)= integral(@(x) nBlock(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(2,1)= integral(@(x) pBlock(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Left 2x4 block
matrixDelN(3,1) = integral(@(x) nBlock(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(4,1) = integral(@(x) nBlock(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelN(5,1) = integral(@(x) pBlock(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(6,1) = integral(@(x) pBlock(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Right 4x2 block
matrixDelN(1,3) = integral(@(x) nBlock(x).*f2((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(1,4) = integral(@(x) nBlock(x).*f4((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelN(2,3) = integral(@(x) pBlock(x).*f2((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(2,4) = integral(@(x) pBlock(x).*f4((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelN(phi(n,r),xx)
matrixDelN(3,3) = matrixDelN(3,3) + integral(@(x) nBlock(x).*f2((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(4,3) = matrixDelN(4,3) + integral(@(x) nBlock(x).*f2((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelN(psi(n,r),xx)
matrixDelN(3,4) = matrixDelN(3,4) + integral(@(x) nBlock(x).*f4((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(4,4) = matrixDelN(4,4) + integral(@(x) nBlock(x).*f4((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelN(phi(p,r),xx)
matrixDelN(5,3) = matrixDelN(5,3) + integral(@(x) pBlock(x).*f2((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(6,3) = matrixDelN(6,3) + integral(@(x) pBlock(x).*f2((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelN(psi(p,r),xx)
matrixDelN(5,4) = matrixDelN(5,4) + integral(@(x) pBlock(x).*f4((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(6,4) = matrixDelN(6,4) + integral(@(x) pBlock(x).*f4((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

    k=numIntervals;
    nBlock = @(x) 0;%Cn3.*((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))) - ((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2)-1/taun;
    pBlock = @(x) -Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;
           
%Handle 1st interval 'half case'
%Bottom right 2x2 block
matrixDelN(last,last-1) =   integral(@(x) pBlock(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-1,last-1) = integral(@(x) nBlock(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Bottom 2x4 block
matrixDelN(last,last-4) = integral(@(x) pBlock(x).*f3((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last,last-5) = integral(@(x) pBlock(x).*f1((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelN(last-1,last-4) = integral(@(x) nBlock(x).*f3((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-1,last-5) = integral(@(x) nBlock(x).*f1((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Left 4x2 Block
matrixDelN(last-2,last-1) = integral(@(x) pBlock(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-3,last-1) = integral(@(x) pBlock(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelN(last-4,last-1) = integral(@(x) nBlock(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-5,last-1) = integral(@(x) nBlock(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelN(phi(n,r),xx)
matrixDelN(last-5,last-5) = matrixDelN(last-5,last-5) + integral(@(x) nBlock(x).*f1(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-4,last-5) = matrixDelN(last-4,last-5) + integral(@(x) nBlock(x).*f1(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelN(psi(n,r),xx)
matrixDelN(last-5,last-4) = matrixDelN(last-5,last-4) + integral(@(x) nBlock(x).*f3(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-4,last-4) = matrixDelN(last-4,last-4) + integral(@(x) nBlock(x).*f3(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelN(phi(p,r),xx)
matrixDelN(last-3,last-5) = matrixDelN(last-3,last-5) + integral(@(x) pBlock(x).*f1(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-2,last-5) = matrixDelN(last-2,last-5) + integral(@(x) pBlock(x).*f1(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelN(psi(p,r),xx)
matrixDelN(last-3,last-4) = matrixDelN(last-3,last-4) + integral(@(x) pBlock(x).*f3(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelN(last-2,last-4) = matrixDelN(last-2,last-4) + integral(@(x) pBlock(x).*f3(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
  
Matrix3Finaln = matrixDelN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of n(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with delp(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

matrixDelP = sparse(4*(numIntervals),4*(numIntervals));

for k = 1:numIntervals-2
   k=k+1;
    %2*Cn3*(pk')^2*nk/(pk)^3-Cn3*pk'' nk/(pk)^2- Cn2*pk' nk'/(pk)^2 
    nBlock = @(x) -Cn2.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2;
    %Cp3*((nk'/nk)^2-nk''/nk)-1/taup;
    pBlock = @(x) 0;%Cp3*(((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))-1/taup;
    k=k-1;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 1  - (Upper left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelP(phi(n,r),xx) - nEq
	matrixDelP(4*k-1,4*k+1) = matrixDelP(4*k-1,4*k+1) + integral(@(x) nBlock(x).*f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k,4*k+1)   = matrixDelP(4*k,4*k+1) +   integral(@(x) nBlock(x).*f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(n,r),xx) - nEq
    matrixDelP(4*k-1,4*k+2) = matrixDelP(4*k-1,4*k+2) + integral(@(x) nBlock(x).*f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k,4*k+2)   = matrixDelP(4*k,4*k+2) +   integral(@(x) nBlock(x).*f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(phi(p,r),xx) - pEq
	matrixDelP(4*k+1,4*k+1) = matrixDelP(4*k+1,4*k+1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+2,4*k+1) = matrixDelP(4*k+2,4*k+1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(p,r),xx) - pEq
	matrixDelP(4*k+1,4*k+2) = matrixDelP(4*k+1,4*k+2) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+2,4*k+2) = matrixDelP(4*k+2,4*k+2) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 2  - (Upper right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelP(phi(n,r),xx) - nEq
	matrixDelP(4*k-1,4*k+5) = matrixDelP(4*k-1,4*k+5) + integral(@(x) nBlock(x).*f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k,4*k+5)   = matrixDelP(4*k,4*k+5) +   integral(@(x) nBlock(x).*f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(n,r),xx) - nEq
    matrixDelP(4*k-1,4*k+6) = matrixDelP(4*k-1,4*k+6) + integral(@(x) nBlock(x).*f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k,4*k+6)   = matrixDelP(4*k,4*k+6) +   integral(@(x) nBlock(x).*f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(phi(p,r),xx) - pEq
	matrixDelP(4*k+1,4*k+5) = matrixDelP(4*k+1,4*k+5) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+2,4*k+5) = matrixDelP(4*k+2,4*k+5) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(p,r),xx) - pEq
	matrixDelP(4*k+1,4*k+6) = matrixDelP(4*k+1,4*k+6) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+2,4*k+6) = matrixDelP(4*k+2,4*k+6) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 3  - (Bottom left 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelP(phi(n,r),xx) - nEq
	matrixDelP(4*k+3,4*k+1) = matrixDelP(4*k+3,4*k+1) + integral(@(x) nBlock(x).*f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+4,4*k+1) = matrixDelP(4*k+4,4*k+1) + integral(@(x) nBlock(x).*f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(n,r),xx) - nEq
    matrixDelP(4*k+3,4*k+2) = matrixDelP(4*k+3,4*k+2) + integral(@(x) nBlock(x).*f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+4,4*k+2) = matrixDelP(4*k+4,4*k+2) + integral(@(x) nBlock(x).*f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(phi(p,r),xx) - pEq
	matrixDelP(4*k+5,4*k+1) = matrixDelP(4*k+5,4*k+1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+6,4*k+1) = matrixDelP(4*k+6,4*k+1) + integral(@(x) pBlock(x).*f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(p,r),xx) - pEq
	matrixDelP(4*k+5,4*k+2) = matrixDelP(4*k+5,4*k+2) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+6,4*k+2) = matrixDelP(4*k+6,4*k+2) + integral(@(x) pBlock(x).*f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Block 4  - (Bottom right 4x4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %matrixDelP(phi(n,r),xx) - nEq
	matrixDelP(4*k+3,4*k+5) = matrixDelP(4*k+3,4*k+5) + integral(@(x) nBlock(x).*f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+4,4*k+5) = matrixDelP(4*k+4,4*k+5) + integral(@(x) nBlock(x).*f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(n,r),xx) - nEq
    matrixDelP(4*k+3,4*k+6) = matrixDelP(4*k+3,4*k+6) + integral(@(x) nBlock(x).*f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+4,4*k+6) = matrixDelP(4*k+4,4*k+6) + integral(@(x) nBlock(x).*f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(phi(p,r),xx) - pEq
	matrixDelP(4*k+5,4*k+5) = matrixDelP(4*k+5,4*k+5) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+6,4*k+5) = matrixDelP(4*k+6,4*k+5) + integral(@(x) pBlock(x).*f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    %matrixDelP(psi(p,r),xx) - pEq
	matrixDelP(4*k+5,4*k+6) = matrixDelP(4*k+5,4*k+6) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	matrixDelP(4*k+6,4*k+6) = matrixDelP(4*k+6,4*k+6) + integral(@(x) pBlock(x).*f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
   
end    
    k=1;
    nBlock = @(x) -Cn2.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2;
    pBlock = @(x) 0;%Cp3*(((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))-1/taup;
    
%Handle 1st interval 'half case'
%Top left 2x2 block
matrixDelP(1,2)= integral(@(x) nBlock(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(2,2)= integral(@(x) pBlock(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Left 2x4 block
matrixDelP(3,2) = integral(@(x) nBlock(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(4,2) = integral(@(x) nBlock(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelP(5,2) = integral(@(x) pBlock(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(6,2) = integral(@(x) pBlock(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Right 4x2 block
matrixDelP(1,5) = integral(@(x) nBlock(x).*f2((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(1,6) = integral(@(x) nBlock(x).*f4((x-node(1))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrixDelP(2,5) = integral(@(x) pBlock(x).*f2((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(2,6) = integral(@(x) pBlock(x).*f4((x-node(1))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelP(phi(n,r),xx)
matrixDelP(3,5) = matrixDelP(3,5) + integral(@(x) nBlock(x).*f2((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(4,5) = matrixDelP(4,5) + integral(@(x) nBlock(x).*f2((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelP(psi(n,r),xx)
matrixDelP(3,6) = matrixDelP(3,6) + integral(@(x) nBlock(x).*f4((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(4,6) = matrixDelP(4,6) + integral(@(x) nBlock(x).*f4((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelP(phi(p,r),xx)
matrixDelP(5,5) = matrixDelP(5,5) + integral(@(x) pBlock(x).*f2((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(6,5) = matrixDelP(6,5) + integral(@(x) pBlock(x).*f2((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelP(psi(p,r),xx)
matrixDelP(5,6) = matrixDelP(5,6) + integral(@(x) pBlock(x).*f4((x-node(1))).*f2((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(6,6) = matrixDelP(6,6) + integral(@(x) pBlock(x).*f4((x-node(1))).*f4((x-node(1))),node(1),node(2),'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

    k=numIntervals;
    nBlock = @(x) -Cn2.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2;
    pBlock = @(x) 0;%Cp3*(((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))-1/taup;
    
%Handle 1st interval 'half case'
%Bottom right 2x2 block
matrixDelP(last,last) =   integral(@(x) pBlock(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-1,last) = integral(@(x) nBlock(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Bottom 2x4 block
matrixDelP(last,last-2) = integral(@(x) pBlock(x).*f3((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last,last-3) = integral(@(x) pBlock(x).*f1((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelP(last-1,last-2) = integral(@(x) nBlock(x).*f3((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-1,last-3) = integral(@(x) nBlock(x).*f1((x-node(end)+dx)).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Left 4x2 Block
matrixDelP(last-2,last) = integral(@(x) pBlock(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-3,last) = integral(@(x) pBlock(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrixDelP(last-4,last) = integral(@(x) nBlock(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-5,last) = integral(@(x) nBlock(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%Overlapping 4x4 block
%matrixDelP(phi(n,r),xx)
matrixDelP(last-5,last-3) = matrixDelP(last-5,last-3) + integral(@(x) nBlock(x).*f1(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-4,last-3) = matrixDelP(last-4,last-3) + integral(@(x) nBlock(x).*f1(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelP(psi(n,r),xx)
matrixDelP(last-5,last-2) = matrixDelP(last-5,last-2) + integral(@(x) nBlock(x).*f3(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-4,last-2) = matrixDelP(last-4,last-2) + integral(@(x) nBlock(x).*f3(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelP(phi(p,r),xx)
matrixDelP(last-3,last-3) = matrixDelP(last-3,last-3) + integral(@(x) pBlock(x).*f1(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-2,last-3) = matrixDelP(last-2,last-3) + integral(@(x) pBlock(x).*f1(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%matrixDelP(psi(p,r),xx)
matrixDelP(last-3,last-2) = matrixDelP(last-3,last-2) + integral(@(x) pBlock(x).*f3(x-node(end)+dx).*f1(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrixDelP(last-2,last-2) = matrixDelP(last-2,last-2) + integral(@(x) pBlock(x).*f3(x-node(end)+dx).*f3(x-node(end)+dx),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
   
Matrix3Finalp = matrixDelP;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of delp(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Residual terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Residual from n'' term
%Should only have values for Cns and Dns, not Ans and Bns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Term is Cn1*n''(x)

%Generate matrix of test functions against RHS of equation
Rn1 = zeros(4*(numIntervals),1);

for k = 1:numIntervals-2
    k=k+1;
    nBlock =@(x) Cn1*(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)));
   k=k-1;
   
    Rn1(4*k-1) =      Rn1(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn1(4*k) =        Rn1(4*k)   +        integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    Rn1(4*k+3) =      Rn1(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn1(4*k+4) =      Rn1(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end
k=1;
nBlock =@(x) Cn1*(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)));
  
%1st interval edge cases
Rn1(1) = Rn1(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rn1(3) = Rn1(3) + integral(@(x) f2((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rn1(4) = Rn1(4) + integral(@(x) f4((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;
k=numIntervals;
 nBlock =@(x) Cn1*(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)));
   
%last interval edge cases
Rn1(last-1) = Rn1(last-1)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rn1(last-5) = Rn1(last-5) + integral(@(x) f1((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rn1(last-4) = Rn1(last-4) + integral(@(x) f3((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%%%%%%%%%%%%%%%%%%%%%
%End residual n'' term
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
%Begin residual p'' term
%%%%%%%%%%%%%%%%%%%%%%
%Term is Cp1*p''

Rp1 = zeros(4*(numIntervals),1);


for k = 1:numIntervals-2
    k=k+1;
    pBlock =@(x) Cp1*(Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)));
    k=k-1;
    %pGuess terms
	Rp1(4*k+1) =    Rp1(4*k+1) +      integral(@(x) f1((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rp1(4*k+2) =    Rp1(4*k+2) +      integral(@(x) f3((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    Rp1(4*k+5) =    Rp1(4*k+5) +      integral(@(x) f2((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rp1(4*k+6) =    Rp1(4*k+6) +      integral(@(x) f4((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
       
end

k=1; 
pBlock =@(x) Cp1*(Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)));
    
%1st interval edge cases
Rp1(2) = Rp1(2) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rp1(5) = Rp1(5) + integral(@(x) f2((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rp1(6) = Rp1(6) + integral(@(x) f4((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

k=numIntervals; 
pBlock =@(x) Cp1*(Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)));

%Last interval edge cases
Rp1(last)   = Rp1(last)   + integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rp1(last-3) = Rp1(last-3) + integral(@(x) f1((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rp1(last-2) = Rp1(last-2) + integral(@(x) f3((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End residual p'' term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Residual from n' term
%Should only have values for Cns and Dns, not Ans and Bns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Term is -Cn2*pk'*nk'/pk

%Generate matrix of test functions against RHS of equation
Rn2 = zeros(4*(numIntervals),1);

for k = 1:numIntervals-2
    k=k+1;
    nBlock =@(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    k=k-1;
    
    Rn2(4*k-1) =      Rn2(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn2(4*k) =        Rn2(4*k)   +        integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    Rn2(4*k+3) =      Rn2(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn2(4*k+4) =      Rn2(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end

k=1;
nBlock =@(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    

%1st interval edge cases
Rn2(1) = Rn2(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rn2(3) = Rn2(3) + integral(@(x) f2((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rn2(4) = Rn2(4) + integral(@(x) f4((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

k=numIntervals;
nBlock =@(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));

%last interval edge cases
Rn2(last-1) = Rn2(last-1)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rn2(last-5) = Rn2(last-5) + integral(@(x) f1((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rn2(last-4) = Rn2(last-4) + integral(@(x) f3((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%%%%%%%%%%%%%%%%%%%%%
%End residual n' term
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
%Begin residual p' term
%%%%%%%%%%%%%%%%%%%%%%
%Term is -Cp2*n'*p'/n

Rp2 = zeros(4*(numIntervals),1);

for k = 1:numIntervals-2
    k=k+1;
    pBlock =@(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
    k=k-1;
    %pGuess terms
	Rp2(4*k+1) =    Rp2(4*k+1) +      integral(@(x) f1((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rp2(4*k+2) =    Rp2(4*k+2) +      integral(@(x) f3((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    Rp2(4*k+5) =    Rp2(4*k+5) +      integral(@(x) f2((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rp2(4*k+6) =    Rp2(4*k+6) +      integral(@(x) f4((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
       
end

 k=1;
pBlock =@(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
 
%1st interval edge cases
Rp2(2) = Rp2(2) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rp2(5) = Rp2(5) + integral(@(x) f2((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rp2(6) = Rp2(6) + integral(@(x) f4((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

 k=numIntervals;
pBlock =@(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));

%Last interval edge cases
Rp2(last)   = Rp2(last)   + integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rp2(last-3) = Rp2(last-3) + integral(@(x) f1((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rp2(last-2) = Rp2(last-2) + integral(@(x) f3((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End residual p' term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Residual from n term
%Should only have values for Cns and Dns, not Ans and Bns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Term is (Cn3*(pk''/pk - (pk'/pk)^2)-1/taun)*nk

%Generate matrix of test functions against RHS of equation
Rn3 = zeros(4*(numIntervals),1);
% 
% for k = 1:numIntervals-2
%     k=k+1;
%     nBlock =@(x) (Cn3*(((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))) - (((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))))).^2)-1/taun).*((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
%     k=k-1;
%     
%     Rn3(4*k-1) =      Rn3(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn3(4*k) =        Rn3(4*k)   +        integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%     Rn3(4*k+3) =      Rn3(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn3(4*k+4) =      Rn3(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% end
% 
% k=1;
% nBlock =@(x) (Cn3*(((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))) - (((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))))).^2)-1/taun).*((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
%     
% %1st interval edge cases
% Rn3(1) = Rn3(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rn3(3) = Rn3(3) + integral(@(x) f2((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rn3(4) = Rn3(4) + integral(@(x) f4((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% 
% %Handle last interval 'half case'
% last = 4*numIntervals;
% 
% k=numIntervals;
% nBlock =@(x) (Cn3*(((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))) - (((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))))).^2)-1/taun).*((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
%  
% %last interval edge cases
% Rn3(last-1) = Rn3(last-1)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rn3(last-5) = Rn3(last-5) + integral(@(x) f1((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rn3(last-4) = Rn3(last-4) + integral(@(x) f3((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%%%%%%%%%%%%%%%%%%%%%
%End residual n term
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
%Begin residual p term
%%%%%%%%%%%%%%%%%%%%%%
%Term is (Cp3*((nk'/nk)^2-nk''/nk)-1/taup)*pk

Rp3 = zeros(4*(numIntervals),1);
% 
% for k = 1:numIntervals-2
%     k=k+1;
%     pBlock =@(x) (Cp3*((((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))).^2 - ((Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))))-1/taup).*((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))));
%     k=k-1;
%     
%     %pGuess terms
% 	Rp3(4*k+1) =    Rp3(4*k+1) +      integral(@(x) f1((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rp3(4*k+2) =    Rp3(4*k+2) +      integral(@(x) f3((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%     Rp3(4*k+5) =    Rp3(4*k+5) +      integral(@(x) f2((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rp3(4*k+6) =    Rp3(4*k+6) +      integral(@(x) f4((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%        
% end
% 
% k=1;
% pBlock =@(x) (Cp3*((((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))).^2 - ((Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))))-1/taup).*((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))));
%     
% %1st interval edge cases
% Rp3(2) = Rp3(2) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rp3(5) = Rp3(5) + integral(@(x) f2((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rp3(6) = Rp3(6) + integral(@(x) f4((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% 
% %Handle last interval 'half case'
% last = 4*numIntervals;
% 
% k=numIntervals;
% pBlock =@(x) (Cp3*((((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))).^2 - ((Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))))-1/taup).*((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))));
% 
% %Last interval edge cases
% Rp3(last)   = Rp3(last)   + integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rp3(last-3) = Rp3(last-3) + integral(@(x) f1((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rp3(last-2) = Rp3(last-2) + integral(@(x) f3((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End residual p term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jacobian = Matrix1Finaln + Matrix1Finalp + Matrix2Finaln + Matrix2Finalp + Matrix3Finaln + Matrix3Finalp;
output = RHSOut-(-Rn1-Rn2-Rn3-Rp1-Rp2-Rp3);
end