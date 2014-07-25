function output = callFunction(inputCoeffs)

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
p1=-2.19258*10^11;
p2=8.16038*10^13;
p3=-8*10^13;
gammap1=100000;
gammap2=100;
gammap3=1;

%Symbolic form of pGuess
pGuessOrig = @(x) p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(-gammap3*x);

n1=-1.79577*10^10;
n2=1.29167*10^11;
n3=10^12;
gamman1=100000;
gamman2=10000;
gamman3=1;

%Symbolic form of nGuess
nGuessOrig = @(x) n1*exp(-gamman1*x)+n2*exp(-gamman2*x)+n3*exp(-gamman3*x);

nprimeOrig = @(x) -(n1*gamman1*exp(-gamman1*x) + n2*gamman2*exp(-gamman2*x) + n3*gamman3*exp(-gamman3*x));
pprimeOrig = @(x) -(p1*gammap1*exp(-gammap1*x) + p2*gammap2*exp(-gammap2*x) + p3*gammap3*exp(-gammap3*x));

nprimeprimeOrig = @(x) (n1*gamman1*gamman1*exp(-gamman1*x) + n2*gamman2*gamman2*exp(-gamman2*x) + n3*gamman3*gamman3*exp(-gamman3*x));
pprimeprimeOrig = @(x) (p1*gammap1*gammap1*exp(-gammap1*x) + p2*gammap2*gammap2*exp(-gammap2*x) + p3*gammap3*gammap3*exp(-gammap3*x));

cn1 = Dn - ksiN*T*mun/q;
cn2 = ksiN*T*mun/q;

cp1 = ksiP*T*mup/q - Dp;
cp2 = -ksiP*T*mup/q;


Gopn = @(x) cn1.*(n1*(gamman1*gamman1)*exp(-gamman1*x) + n2*(gamman2*gamman2)*exp(-gamman2*x) + n3*(gamman3*gamman3)*exp(-gamman3*x))    + cn2.*pprimeOrig(x)./pGuessOrig(x).*(-n1*(gamman1)*exp(-gamman1*x) - n2*(gamman2)*exp(-gamman2*x) - n3*(gamman3)*exp(-gamman3*x));%+(ksiN*T*mun/q*(pprimeprimeOrig(x)./pGuessOrig(x) - (pprimeOrig(x)./pGuessOrig(x)).^2) - 1/taun).*(n1*exp(-gamman1*x) + n2*exp(-gamman2*x) + n3*exp(-gamman3*x));     
Gopp = @(x) cp1.*(p1*(gammap1*gammap1)*exp(-gammap1*x) + p2*(gammap2*gammap2)*exp(-gammap2*x) + p3*(gammap3*gammap3)*exp(-gammap3*x))    + cp2.*nprimeOrig(x)./nGuessOrig(x).*(-p1*(gammap1)*exp(-gammap1*x) - p2*(gammap2)*exp(-gammap2*x) - p3*(gammap3)*exp(-gammap3*x));%+(ksiP*T*mup/q*((nprimeOrig(x)./nGuessOrig(x)).^2 - nprimeprimeOrig(x)./nGuessOrig(x)) - 1/taup).*(p1*exp(-gammap1*x) + p2*exp(-gammap2*x) + p3*exp(-gammap3*x));     
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
%BCs, copied from driver function, not calculated because it would be a
%PITA
%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphan=454.548;
alphap=10000;
betan=-454.546;
betap=-10000;

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
% Rn3right = zeros(4*(numIntervals),1);
% Rn3left = zeros(4*(numIntervals),1);
% for k = 1:numIntervals-2
%     k=k+1;
%     nBlock =@(x) (Cn3*(((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))) - (((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))))).^2)-1/taun).*((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
%     k=k-1;
%     
%     Rn3(4*k-1) =      Rn3(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn3(4*k) =        Rn3(4*k)   +        integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%     Rn3(4*k+3) =      Rn3(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn3(4*k+4) =      Rn3(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%     
%     Rn3right(4*k-1) =    integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn3right(4*k) =      integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%     Rn3left(4*k+3) =     integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn3left(4*k+4) =     integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
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
% 
% 
% %%%%%%%%%%%%%%%%%%%%%
% %End residual n term
% %%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% %Begin residual p term
% %%%%%%%%%%%%%%%%%%%%%%
% %Term is (Cp3*((nk'/nk)^2-nk''/nk)-1/taup)*pk
% 
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

output = RHSOut-Rn1-Rn2-Rn3-Rp1-Rp2-Rp3;
end