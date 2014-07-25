clear all
close all
t=cputime;
%profile on


%Initial Newton method calculation utilizing a known test case
%Only n'',n', p'', p' terms
%No noise


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
taun = 10^-12; %(*electron lifetime in s*)
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
rand=0.95;

AbsError=10^-18;
RelError=10^-6;

filepath=strcat('newton v9/rand',num2str(rand),'/',num2str(numIntervals),'terms');
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


Gopn = @(x) cn1.*(n1*(gamman1*gamman1)*exp(-gamman1*x) + n2*(gamman2*gamman2)*exp(-gamman2*x) + n3*(gamman3*gamman3)*exp(-gamman3*x))    +(ksiN*T*mun/q*(pprimeprimeOrig(x)./pGuessOrig(x) - (pprimeOrig(x)./pGuessOrig(x)).^2) - 1/taun).*(n1*exp(-gamman1*x) + n2*exp(-gamman2*x) + n3*exp(-gamman3*x));%+cn2.*pprime(x)./pGuess(x).*(-n1*(gamman1)*exp(-gamman1*x) - n2*(gamman2)*exp(-gamman2*x) - n3*(gamman3)*exp(-gamman3*x))     
Gopp = @(x) cp1.*(p1*(gammap1*gammap1)*exp(-gammap1*x) + p2*(gammap2*gammap2)*exp(-gammap2*x) + p3*(gammap3*gammap3)*exp(-gammap3*x))    +(ksiP*T*mup/q*((nprimeOrig(x)./nGuessOrig(x)).^2 - nprimeprimeOrig(x)./nGuessOrig(x)) - 1/taup).*(p1*exp(-gammap1*x) + p2*exp(-gammap2*x) + p3*exp(-gammap3*x));%+cp2.*nprime(x)./nGuess(x).*(-p1*(gammap1)*exp(-gammap1*x) - p2*(gammap2)*exp(-gammap2*x) - p3*(gammap3)*exp(-gammap3*x))     

%Adjust pGuess to include error
%Adjust pGuess to include error

%'Noisified' terms
%n3*rand,gamman1,gamman2,gamman3 set;n1 and n2 calculated from those
%numbers using BCs
%Calcs done in Mathematica

n1rand=-1.70787*10^10;
n2rand=1.22709*10^11;
n3rand=10^12*(rand);

nGuess = @(x) (n1rand*exp(-gamman1*x)+n2rand*exp(-gamman2*x)+n3rand*exp(-gamman3*x));
nprime = @(x) -(n1rand*gamman1*exp(-gamman1*x) + n2rand*gamman2*exp(-gamman2*x) + n3rand*gamman3*exp(-gamman3*x));
nprimeprime = @(x) n1rand*gamman1*gamman1*exp(-gamman1*x) + n2rand*gamman2*gamman2*exp(-gamman2*x) + n3rand*gamman3*gamman3*exp(-gamman3*x);

%p3*rand,gammap1,gammap2,gammap3 set;p1 and p2 calculated from those
%numbers using BCs
%Calcs done in Mathematica

p1rand=-2.08295*10^11;
p2rand=7.75236*10^13;
p3rand=-8*10^13*rand;

pGuess = @(x) (p1rand*exp(-gammap1*x)+p2rand*exp(-gammap2*x)+p3rand*exp(-gammap3*x));
pprime = @(x) -(p1rand*gammap1*exp(-gammap1*x) + p2rand*gammap2*exp(-gammap2*x) + p3rand*gammap3*exp(-gammap3*x));
pprimeprime =@(x) p1rand*gammap1*gammap1*exp(-gammap1*x) + p2rand*gammap2*gammap2*exp(-gammap2*x) + p3rand*gammap3*gammap3*exp(-gammap3*x);

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
alphan = nprime(0)/nGuess(0);

%beta is boundary condition term, n'(L) = beta n(L)
betan= nprime(totalLength)/nGuess(totalLength);

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alphap = pprime(0)/pGuess(0);

%beta is boundary condition term, n'(L) = beta n(L)
betap= pprime(totalLength)/pGuess(totalLength);


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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Add in BC values omitted from calcs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cnjs(1,1) = Coeffsn(1);
% Dnjs(1,1) = alpha*Coeffsn(1);
% Cnjs(1,numIntervals+1) = Coeffsn(length(Coeffsn));
% Dnjs(1,numIntervals+1) = beta*Coeffsn(length(Coeffsn));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Add in BC values omitted from calcs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cpjs(1,1) = Coeffsp(1);
% Dpjs(1,1) = alphap*Coeffsp(1);
% Cpjs(1,numIntervals+1) = Coeffsp(length(Coeffsp));
% Dpjs(1,numIntervals+1) = betap*Coeffsp(length(Coeffsp));

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
% 
% fig4 = figure;
% plot(indexn,Dnjs,'x',node,nprime(node),indexn,nprimeOrig(indexn)),title('nprime (recomposition)'),legend('Dnjs','Error','Orginial');
% saveas(fig4,strcat(filepath,'/nprime_recomposition.jpg'))

fig5 = figure;
plot(indexn,Cpjs,'x',node,pGuess(node),indexn,pOrig(indexn)),title('pGuess (recomposition)'),legend('Cpjs','Error','Orginial');
saveas(fig5,strcat(filepath,'/pGuess_recomposition.jpg'))
% 
% fig5 = figure;
% plot(indexn,Dpjs,'x',node,pprime(node),indexn,pprimeOrig(indexn)),title('pprime (recomposition)'),legend('Dpjs','Error','Orginial');
% saveas(fig5,strcat(filepath,'/pprime_recomposition.jpg'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General BCs for calculated solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alphan = nprimeOrig(0)/nOrig(0);%Sr/Dn;

%beta is boundary condition term, n'(L) = beta n(L)
betan = nprimeOrig(totalLength)/nOrig(totalLength);%-Sr/Dn;

%alpha is boundary conditions term, p'(0)=alpha*p(0)
alphap = pprimeOrig(0)/pOrig(0);%Sr/Dp;
%beta is boundary condition term, p'(L) = beta p(L)
betap = pprimeOrig(totalLength)/pOrig(totalLength);%-Sr/Dp;

%pGuess = matlabFunction(pGuess_syms);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End BCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cn1=Dn - ksiN*T*mun/q;
Cn2=ksiN*T*mun/q;
Cn3=ksiN*T*mun/q;
Cp1=ksiP*T*mup/q - Dp;
Cp2=-ksiP*T*mup/q;
Cp3=ksiP*T*mup/q;

for iterationCount=1:maxIterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates for RHS of equation (forcing function times test functions)
%This calculation is performed outside the iterative loop, as it will not
%change over the course of changing p(x), n(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate matrix of test functions against RHS of equation
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End RHS calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with (deln)''(x) term
% %Performed outside iterative loop because this term is not influenced by
% %p(x) or n(x).  Should not change over calculations.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Front end coefficient of n''(x) term

matrixDelNpp = sparse(4*(numIntervals),4*(numIntervals));

for k = 1:numIntervals-2
   k=k+1;
    %Cn1
    nBlock = @(x) Cn1;
    %-Cp3*pk/nk  
    pBlock = @(x) -Cp3.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
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
    pBlock = @(x) -Cp3.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
      
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
    pBlock = @(x) -Cp3.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
    
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
    nBlock = @(x) Cn3.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
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
    nBlock = @(x) Cn3.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
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
    nBlock = @(x) Cn3.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
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
    nBlock = @(x) 0;%Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    %Cp2 p'k/nk + 2*Cn3*(n'k) pk/(nk)^2   
    pBlock = @(x) 2.*Cp3.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;%Cp2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
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
    nBlock = @(x) 0;%Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) 2.*Cp3.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;%Cp2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
       
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
    nBlock = @(x) 0;%Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) 2.*Cp3.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;%Cp2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
      
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
    nBlock = @(x) -2.*Cn3.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2);%Cn2.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))/(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    %Cp2*nk'/nk   
    pBlock = @(x) 0;%Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
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
    nBlock = @(x) -2.*Cn3.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2);%Cn2.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))/(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) 0;%Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
         
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
    nBlock = @(x) -2.*Cn3.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2);%Cn2.*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))/(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    pBlock = @(x) 0;%Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
        
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
    nBlock = @(x) Cn3.*((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))) - ((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2)-1/taun;
    %Cp3*nk''pk/nk^2-Cp2*nk' pk'/nk^2 - 2Cp3(nk')^2*pk/(nk)^3
    pBlock = @(x) Cp3*(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))).*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2 - 2*Cp3*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).^2.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^3;%- Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;
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
    nBlock = @(x) Cn3.*((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))) - ((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2)-1/taun;
    pBlock = @(x) Cp3*(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))).*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2 - 2*Cp3*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).^2.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^3;%- Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;
    
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
    nBlock = @(x) Cn3.*((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))) - ((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2)-1/taun;
    pBlock = @(x) Cp3*(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))).*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2 - 2*Cp3*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).^2.*(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^3;%- Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))).^2;
           
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
    nBlock = @(x) 2.*Cn3.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).^2.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^3)-Cn3.*(Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))).*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2);%-Cn2.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2; 
    %Cp3*((nk'/nk)^2-nk''/nk)-1/taup;
    pBlock = @(x) Cp3*(((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))-1/taup;
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
    nBlock = @(x) 2.*Cn3.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).^2.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^3)-Cn3.*(Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))).*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2);%-Cn2.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2; 
    pBlock = @(x) Cp3*(((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))-1/taup;
    
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
    nBlock = @(x) 2.*Cn3.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).^2.*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^3)-Cn3.*(Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))).*(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2);%-Cn2.*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))).^2; 
    pBlock = @(x) Cp3*(((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))-1/taup;
    
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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Residual from n' term
% %Should only have values for Cns and Dns, not Ans and Bns
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Term is -Cn2*pk'*nk'/pk
% 
% %Generate matrix of test functions against RHS of equation
% Rn2 = zeros(4*(numIntervals),1);
% 
% for k = 1:numIntervals-2
%     k=k+1;
%     nBlock =@(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
%     k=k-1;
%     
%     Rn2(4*k-1) =      Rn2(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn2(4*k) =        Rn2(4*k)   +        integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%     Rn2(4*k+3) =      Rn2(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rn2(4*k+4) =      Rn2(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% end
% 
% k=1;
% nBlock =@(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
%     
% 
% %1st interval edge cases
% Rn2(1) = Rn2(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rn2(3) = Rn2(3) + integral(@(x) f2((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rn2(4) = Rn2(4) + integral(@(x) f4((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% 
% %Handle last interval 'half case'
% last = 4*numIntervals;
% 
% k=numIntervals;
% nBlock =@(x) Cn2*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))).*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
% 
% %last interval edge cases
% Rn2(last-1) = Rn2(last-1)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rn2(last-5) = Rn2(last-5) + integral(@(x) f1((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rn2(last-4) = Rn2(last-4) + integral(@(x) f3((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% 
% %%%%%%%%%%%%%%%%%%%%%
% %End residual n' term
% %%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%
% %Begin residual p' term
% %%%%%%%%%%%%%%%%%%%%%%
% %Term is -Cp2*n'*p'/n
% 
% Rp2 = zeros(4*(numIntervals),1);
% 
% for k = 1:numIntervals-2
%     k=k+1;
%     pBlock =@(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
%     k=k-1;
%     %pGuess terms
% 	Rp2(4*k+1) =    Rp2(4*k+1) +      integral(@(x) f1((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rp2(4*k+2) =    Rp2(4*k+2) +      integral(@(x) f3((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%     Rp2(4*k+5) =    Rp2(4*k+5) +      integral(@(x) f2((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
% 	Rp2(4*k+6) =    Rp2(4*k+6) +      integral(@(x) f4((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
%        
% end
% 
%  k=1;
% pBlock =@(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
%  
% %1st interval edge cases
% Rp2(2) = Rp2(2) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rp2(5) = Rp2(5) + integral(@(x) f2((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% Rp2(6) = Rp2(6) + integral(@(x) f4((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
% 
% %Handle last interval 'half case'
% last = 4*numIntervals;
% 
%  k=numIntervals;
% pBlock =@(x) Cp2*(Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))).*(Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
% 
% %Last interval edge cases
% Rp2(last)   = Rp2(last)   + integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rp2(last-3) = Rp2(last-3) + integral(@(x) f1((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% Rp2(last-2) = Rp2(last-2) + integral(@(x) f3((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End residual p' term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Residual from n term
%Should only have values for Cns and Dns, not Ans and Bns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Term is (Cn3*(pk''/pk - (pk'/pk)^2)-1/taun)*nk

%Generate matrix of test functions against RHS of equation
Rn3 = zeros(4*(numIntervals),1);

Rn3right = zeros(4*(numIntervals),1);
Rn3left = zeros(4*(numIntervals),1);
for k = 1:numIntervals-2
    k=k+1;
    nBlock =@(x) (Cn3*(((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))) - (((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))))).^2)-1/taun).*((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
    k=k-1;
    
    Rn3(4*k-1) =      Rn3(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn3(4*k) =        Rn3(4*k)   +        integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    Rn3(4*k+3) =      Rn3(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn3(4*k+4) =      Rn3(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    Rn3right(4*k-1) =    integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn3right(4*k) =      integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    Rn3left(4*k+3) =     integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rn3left(4*k+4) =     integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end

k=1;
nBlock =@(x) (Cn3*(((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))) - (((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))))).^2)-1/taun).*((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
    
%1st interval edge cases
Rn3(1) = Rn3(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rn3(3) = Rn3(3) + integral(@(x) f2((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rn3(4) = Rn3(4) + integral(@(x) f4((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

k=numIntervals;
nBlock =@(x) (Cn3*(((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))) - (((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k))))./((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))))).^2)-1/taun).*((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
 
%last interval edge cases
Rn3(last-1) = Rn3(last-1)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rn3(last-5) = Rn3(last-5) + integral(@(x) f1((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rn3(last-4) = Rn3(last-4) + integral(@(x) f3((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


trueRn3 = zeros(4*(numIntervals),1);
trueRn3right = zeros(4*(numIntervals),1);
trueRn3left = zeros(4*(numIntervals),1);
internalDiff= zeros(numIntervals,1);
for k = 1:numIntervals-2
    k=k+1;
    nBlock =@(x) (Cn3*((trueCpjs(k).*f1primeprime(x-node(k)) + trueCpjs(k+1).*f2primeprime(x-node(k)) + trueDpjs(k).*f3primeprime(x-node(k)) + trueDpjs(k+1).*f4primeprime(x-node(k)))./(trueCpjs(k).*f1(x-node(k)) + trueCpjs(k+1).*f2(x-node(k)) + trueDpjs(k).*f3(x-node(k)) + trueDpjs(k+1).*f4(x-node(k))) - ((trueCpjs(k).*f1prime(x-node(k)) + trueCpjs(k+1).*f2prime(x-node(k)) + trueDpjs(k).*f3prime(x-node(k)) + trueDpjs(k+1).*f4prime(x-node(k)))./(trueCpjs(k).*f1(x-node(k)) + trueCpjs(k+1).*f2(x-node(k)) + trueDpjs(k).*f3(x-node(k)) + trueDpjs(k+1).*f4(x-node(k)))).^2)-1/taun).*((trueCnjs(k).*f1(x-node(k)) + trueCnjs(k+1).*f2(x-node(k)) + trueDnjs(k).*f3(x-node(k)) + trueDnjs(k+1).*f4(x-node(k))));
    k=k-1;
    
    trueRn3(4*k-1) =      trueRn3(4*k-1) +        integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	trueRn3(4*k) =        trueRn3(4*k)   +        integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    trueRn3(4*k+3) =      trueRn3(4*k+3) +        integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	trueRn3(4*k+4) =      trueRn3(4*k+4) +        integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    trueRn3right(4*k-1) =    integral(@(x) f1((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	trueRn3right(4*k) =      integral(@(x) f3((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    trueRn3left(4*k+3) =     integral(@(x) f2((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	trueRn3left(4*k+4) =     integral(@(x) f4((x-node(k+1))).*nBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    internalDiff(k) = (trueCpjs(k).*f1primeprime(node(k)) + trueCpjs(k+1).*f2primeprime(node(k)) + trueDpjs(k).*f3primeprime(node(k)) + trueDpjs(k+1).*f4primeprime(node(k)))./(trueCpjs(k).*f1(node(k)) + trueCpjs(k+1).*f2(node(k)) + trueDpjs(k).*f3(node(k)) + trueDpjs(k+1).*f4(node(k))) - ((trueCpjs(k).*f1prime(node(k)) + trueCpjs(k+1).*f2prime(node(k)) + trueDpjs(k).*f3prime(node(k)) + trueDpjs(k+1).*f4prime(node(k)))./(trueCpjs(k).*f1(node(k)) + trueCpjs(k+1).*f2(node(k)) + trueDpjs(k).*f3(node(k)) + trueDpjs(k+1).*f4(node(k)))).^2;
end
 internalDiffMod = internalDiff-1/taun;
k=1;
nBlock =@(x) (Cn3*((trueCpjs(k).*f1primeprime(x-node(k)) + trueCpjs(k+1).*f2primeprime(x-node(k)) + trueDpjs(k).*f3primeprime(x-node(k)) + trueDpjs(k+1).*f4primeprime(x-node(k)))./(trueCpjs(k).*f1(x-node(k)) + trueCpjs(k+1).*f2(x-node(k)) + trueDpjs(k).*f3(x-node(k)) + trueDpjs(k+1).*f4(x-node(k))) - ((trueCpjs(k).*f1prime(x-node(k)) + trueCpjs(k+1).*f2prime(x-node(k)) + trueDpjs(k).*f3prime(x-node(k)) + trueDpjs(k+1).*f4prime(x-node(k)))./(trueCpjs(k).*f1(x-node(k)) + trueCpjs(k+1).*f2(x-node(k)) + trueDpjs(k).*f3(x-node(k)) + trueDpjs(k+1).*f4(x-node(k)))).^2)-1/taun).*((trueCnjs(k).*f1(x-node(k)) + trueCnjs(k+1).*f2(x-node(k)) + trueDnjs(k).*f3(x-node(k)) + trueDnjs(k+1).*f4(x-node(k))));
       
%1st interval edge cases
trueRn3(1) = trueRn3(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
trueRn3(3) = trueRn3(3) + integral(@(x) f2((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
trueRn3(4) = trueRn3(4) + integral(@(x) f4((x-node(1))).*nBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

k=numIntervals;
nBlock =@(x) (Cn3*((trueCpjs(k).*f1primeprime(x-node(k)) + trueCpjs(k+1).*f2primeprime(x-node(k)) + trueDpjs(k).*f3primeprime(x-node(k)) + trueDpjs(k+1).*f4primeprime(x-node(k)))./(trueCpjs(k).*f1(x-node(k)) + trueCpjs(k+1).*f2(x-node(k)) + trueDpjs(k).*f3(x-node(k)) + trueDpjs(k+1).*f4(x-node(k))) - ((trueCpjs(k).*f1prime(x-node(k)) + trueCpjs(k+1).*f2prime(x-node(k)) + trueDpjs(k).*f3prime(x-node(k)) + trueDpjs(k+1).*f4prime(x-node(k)))./(trueCpjs(k).*f1(x-node(k)) + trueCpjs(k+1).*f2(x-node(k)) + trueDpjs(k).*f3(x-node(k)) + trueDpjs(k+1).*f4(x-node(k)))).^2)-1/taun).*((trueCnjs(k).*f1(x-node(k)) + trueCnjs(k+1).*f2(x-node(k)) + trueDnjs(k).*f3(x-node(k)) + trueDnjs(k+1).*f4(x-node(k))));
    
%last interval edge cases
trueRn3(last-1) = trueRn3(last-1)+ integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
trueRn3(last-5) = trueRn3(last-5) + integral(@(x) f1((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
trueRn3(last-4) = trueRn3(last-4) + integral(@(x) f3((x-node(end)+dx)).*nBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);


%%%%%%%%%%%%%%%%%%%%%
%End residual n term
%%%%%%%%%%%%%%%%%%%%
% Rn3temp = (Cn3*((pprimeprimeOrig(node)./pOrig(node))-(pprimeOrig(node)./pOrig(node)).^2)-1/taun).*nOrig(node);
% plot(Rn3temp)
% 
% Rp3temp = (Cp3*((nprimeOrig(node)./nOrig(node)).^2-(nprimeprimeOrig(node)./nOrig(node)))-1/taup).*pOrig(node);
% plot(Rp3temp)
%%%%%%%%%%%%%%%%%%%%%%
%Begin residual p term
%%%%%%%%%%%%%%%%%%%%%%
%Term is (Cp3*((nk'/nk)^2-nk''/nk)-1/taup)*pk

Rp3 = zeros(4*(numIntervals),1);

for k = 1:numIntervals-2
    k=k+1;
    pBlock =@(x) (Cp3*((((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))).^2 - ((Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))))-1/taup).*((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))));
    k=k-1;
    
    %pGuess terms
	Rp3(4*k+1) =    Rp3(4*k+1) +      integral(@(x) f1((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rp3(4*k+2) =    Rp3(4*k+2) +      integral(@(x) f3((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    Rp3(4*k+5) =    Rp3(4*k+5) +      integral(@(x) f2((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	Rp3(4*k+6) =    Rp3(4*k+6) +      integral(@(x) f4((x-node(k+1))).*pBlock(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
       
end

k=1;
pBlock =@(x) (Cp3*((((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))).^2 - ((Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))))-1/taup).*((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))));
    
%1st interval edge cases
Rp3(2) = Rp3(2) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rp3(5) = Rp3(5) + integral(@(x) f2((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
Rp3(6) = Rp3(6) + integral(@(x) f4((x-node(1))).*pBlock(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 4*numIntervals;

k=numIntervals;
pBlock =@(x) (Cp3*((((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))))).^2 - ((Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k))))./((Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))))-1/taup).*((Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k))));

%Last interval edge cases
Rp3(last)   = Rp3(last)   + integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rp3(last-3) = Rp3(last-3) + integral(@(x) f1((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
Rp3(last-2) = Rp3(last-2) + integral(@(x) f3((x-node(end)+dx)).*pBlock(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End residual p term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
%Combine Matrix terms, leaving [LHS][coeffs]=[RHS]
%%%%%%%%%%%%%%%%%%%%%%%%

LHSfinal = Matrix1Finaln + Matrix1Finalp + Matrix2Finaln + Matrix2Finalp + Matrix3Finaln + Matrix3Finalp;

RHSFinal = RHSOut-Rn1-Rp1-Rn3-Rp3;%-Rn2-Rp2

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
LHSscaled = D*LHSfinal*D;
RHSscaled = D*RHSFinal;

%Solve simultaneous equations
ScaledCoeffs = LHSscaled\RHSscaled;
Coeffs = D'*ScaledCoeffs;

%Collect Cj's in an array for display purposes
delCnjs = zeros(1,numIntervals+1);
delDnjs = zeros(1,numIntervals+1);
delCpjs = zeros(1,numIntervals+1);
delDpjs = zeros(1,numIntervals+1);

for index=2:length(Cnjs)-1
   delCnjs(1,index) = Coeffs(4*index-5);
   delDnjs(1,index) = Coeffs(4*index-4);
   delCpjs(1,index) = Coeffs(4*index-3);
   delDpjs(1,index) = Coeffs(4*index-2);
end
delCnjs(1) = Coeffs(1);
delDnjs(1) = Coeffs(1)*alphan;
delCpjs(1) = Coeffs(2);
delDpjs(1) = Coeffs(2)*alphap;

delCnjs(numIntervals+1)=Coeffs(length(Coeffs)-1);
delDnjs(numIntervals+1)=Coeffs(length(Coeffs)-1)*betan;
delCpjs(numIntervals+1)=Coeffs(length(Coeffs));
delDpjs(numIntervals+1)=Coeffs(length(Coeffs))*betap;

%% 
%Correct dels
trueDelCnjs = -Cnjs+trueCnjs;
trueDelDnjs = -Dnjs+trueDnjs;
trueDelCpjs = -Cpjs+trueCpjs;
trueDelDpjs = -Dpjs+trueDpjs;

delArray = zeros(4*numIntervals,1);
for index=1:numIntervals-1
    delArray(4*index-1) = trueDelCnjs(index+1);
    delArray(4*index) =   trueDelDnjs(index+1);
    delArray(4*index+1) = trueDelCpjs(index+1);
    delArray(4*index+2) = trueDelDpjs(index+1);
end
delArray(1)=trueDelCnjs(1);
delArray(2)=trueDelCpjs(1);
delArray(end-1)=trueDelCnjs(end);
delArray(end-2)=trueDelCpjs(end);

delArrayStorage{iterationCount}=delArray;

difference = RHSFinal-LHSfinal*delArray;

combo{iterationCount} = LHSfinal*delArray;
RHS{iterationCount} = RHSFinal;

differenceCnjs = [difference(1); difference(3:4:end-2); difference(end-1)];
differenceDnjs = [difference(4:4:end-2);];
differenceCpjs = [difference(2); difference(5:4:end-2); difference(end)];
differenceDpjs = [difference(6:4:end-2)];

differenceStorage{iterationCount} = difference;
differenceCnjsStorage{iterationCount} = differenceCnjs;
differenceDnjsStorage{iterationCount} = differenceDnjs;
differenceCpjsStorage{iterationCount} = differenceCpjs;
differenceDpjsStorage{iterationCount} = differenceDpjs;
%% 
scale=1;
Cnjs = Cnjs+scale*delCnjs;
Cpjs = Cpjs+scale*delCpjs;
Dnjs = Dnjs+scale*delDnjs;
Dpjs = Dpjs+scale*delDpjs;


testCnjs(iterationCount)=norm(delCnjs);
testCpjs(iterationCount)=norm(delCpjs);
testDnjs(iterationCount)=norm(delDnjs);
testDpjs(iterationCount)=norm(delDpjs);

fig4 = figure;
plot(indexn,Cnjs,'x',indexn,nGuess(indexn),indexn,nOrig(indexn));title('n(x) calculated'),legend('Cnjs','nGuess - error','nGuess - original');
saveas(fig4,strcat(filepath,'/n_calculated_iteration',num2str(iterationCount),'.jpg'))

fig6 = figure;
plot(indexn,Dnjs,'x',indexn,nprime(indexn),indexn,nprimeOrig(indexn)),title('nprime calculated'),legend('Dnjs','nprime - error','nprime - original');
saveas(fig6,strcat(filepath,'/nprime_calculated_iteration',num2str(iterationCount),'.jpg'))

% fig5 = figure;
% plot_together(Cnjs,Dnjs,totalLength,numIntervals);title('plot together');
% saveas(fig5,strcat(filepath,'/plot_together_n',num2str(iterationCount),'.jpg'))

fig4 = figure;
plot(indexn,Cpjs,'x',indexn,pGuess(indexn),indexn,pOrig(indexn));title('p(x) calculated'),legend('Cpjs','pGuess - error','pGuess - orignial');
saveas(fig4,strcat(filepath,'/p_calculated_iteration',num2str(iterationCount),'.jpg'))

fig6 = figure;
plot(indexn,Dpjs,'x',indexn,pprime(indexn),indexn,pprimeOrig(indexn)),title('pprime calculated'),legend('Dpjs','pGuess - error','pGuess - orignial');
saveas(fig6,strcat(filepath,'/pprime_calculated_iteration',num2str(iterationCount),'.jpg'))
% 
% fig5 = figure;
% plot_together(Cpjs,Dpjs,totalLength,numIntervals);title('plot together');
% saveas(fig5,strcat(filepath,'/plot_together_p',num2str(iterationCount),'.jpg'))

iterationCount

save(strcat(filepath,'/','iteration',num2str(iterationCount),'.mat'))




end
%Displays runtime of m-file
time = cputime-t



%Store variable values to be used later
save(strcat(filepath,'/',num2str(numIntervals),'terms.mat'))