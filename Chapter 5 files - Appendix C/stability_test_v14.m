clear all
close all
t=cputime;
%profile on

%Adjusts initial pGuess by scaling component exponentials after calculating Gopp, Gopn to test stability

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
dx = intervalLength;
maxIterations=5;
% Store node locations
node = startPoint:dx:totalLength;



%Random noise bounds
randLow=1;%0.99;
randHigh=1;%1.01;

AbsError=10^-16;
RelError=10^-6;

filepath=strcat('stability test v14/rand',num2str(randLow),'/',num2str(numIntervals),'terms');
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
p1=-2.4469*10^11;
p2=1.34584*10^12;
p3=-10^9;
gammap1=100000;
gammap2=10000;
gammap3=1;

%Symbolic form of pGuess
pGuess = @(x) p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(-gammap3*x);

n1=-1.79577*10^10;
n2=1.29167*10^11;
n3=10^12;
gamman1=100000;
gamman2=10000;
gamman3=1;

%Symbolic form of nGuess
nGuess = @(x) n1*exp(-gamman1*x)+n2*exp(-gamman2*x)+n3*exp(-gamman3*x);

nprime = @(x) -(n1*gamman1*exp(-gamman1*x) + n2*gamman2*exp(-gamman2*x) + n3*gamman3*exp(-gamman3*x));
pprime = @(x) -(p1*gammap1*exp(-gammap1*x) + p2*gammap2*exp(-gammap2*x) + p3*gammap3*exp(-gammap3*x));

nprimeprime = @(x) (n1*gamman1*gamman1*exp(-gamman1*x) + n2*gamman2*gamman2*exp(-gamman2*x) + n3*gamman3*gamman3*exp(-gamman3*x));
pprimeprime = @(x) (p1*gammap1*gammap1*exp(-gammap1*x) + p2*gammap2*gammap2*exp(-gammap2*x) + p3*gammap3*gammap3*exp(-gammap3*x));

cn1 = Dn - ksi*T*mun/q;
cn2 = ksi*T*mun/q;

cp1 = ksi*T*mup/q - Dp;
cp2 = -ksi*T*mup/q;


Gopn = @(x) cn1.*(n1*(gamman1*gamman1)*exp(-gamman1*x) + n2*(gamman2*gamman2)*exp(-gamman2*x) + n3*(gamman3*gamman3)*exp(-gamman3*x))    +cn2.*pprime(x)./pGuess(x).*(-n1*(gamman1)*exp(-gamman1*x) - n2*(gamman2)*exp(-gamman2*x) - n3*(gamman3)*exp(-gamman3*x))     +(ksi*T*mun/q*(pprimeprime(x)./pGuess(x) - (pprime(x)./pGuess(x)).^2) - 1/taun).*(n1*exp(-gamman1*x) + n2*exp(-gamman2*x) + n3*exp(-gamman3*x));
Gopp = @(x) cp1.*(p1*(gammap1*gammap1)*exp(-gammap1*x) + p2*(gammap2*gammap2)*exp(-gammap2*x) + p3*(gammap3*gammap3)*exp(-gammap3*x))    +cp2.*nprime(x)./nGuess(x).*(-p1*(gammap1)*exp(-gammap1*x) - p2*(gammap2)*exp(-gammap2*x) - p3*(gammap3)*exp(-gammap3*x))     +(ksi*T*mup/q*((nprime(x)./nGuess(x)).^2 - nprimeprime(x)./nGuess(x)) - 1/taup).*(p1*exp(-gammap1*x) + p2*exp(-gammap2*x) + p3*exp(-gammap3*x));

%Adjust pGuess to include error
pOrig = pGuess;
pprimeOrig = pprime;
nOrig = nGuess;

p1=-2.4469*10^11*(randLow+(randHigh-randLow)*rand(1));
p2=1.34584*10^12*(randLow+(randHigh-randLow)*rand(1));
p3=-10^9*(randLow+(randHigh-randLow)*rand(1));

pGuess = @(x) (p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(-gammap3*x));
pprime = @(x) -(p1*gammap1*exp(-gammap1*x) + p2*gammap2*exp(-gammap2*x) + p3*gammap3*exp(-gammap3*x));

% fig2 = figure;
% plot(node,nGuess(node)), title('nGuess')
% saveas(fig2,strcat(filepath,'/nGuess.jpg'))

%Convert initial guess into form finite element composition
%This will allow it to be usable in later calculations

%Boundary conditions of test function put into form usable for the
%calculations below
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate composition of n(x) in finite element space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alpha is boundary conditions term, n'(0)=alpha*n(0)
alpha = nprime(0)/nGuess(0);

%beta is boundary condition term, n'(L) = beta n(L)
beta= nprime(totalLength)/nGuess(totalLength);

massmatrix = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

massrhs = zeros(2*(numIntervals+1)-2,1);


for k = 1:numIntervals-2
	%%phi_righ*phi_right
	massmatrix(2*k,2*k) = massmatrix(2*k,2*k) + integral(@(x) f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));

	%% phi_right*phi_left
	massmatrix(2*k,2*(k+1)) = massmatrix(2*k,2*(k+1)) + integral(@(x) f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));

       
	%% phi_left*phi_right
	massmatrix(2*(k+1),2*(k)) = massmatrix(2*(k+1),2*(k)) + integral(@(x) f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));

	%% phi_left*phi_left
	massmatrix(2*(k+1),2*(k+1)) = massmatrix(2*(k+1),2*(k+1)) + integral(@(x) f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2)); 

    
    
	%% phi_right*psi_right
	massmatrix(2*k,2*k+1) = massmatrix(2*k,2*k+1) + integral(@(x) f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% phi_right*psi_left
	massmatrix(2*k,2*(k+1)+1) = massmatrix(2*k,2*(k+1)+1) + integral(@(x) f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% phi_left*psi_right
	massmatrix(2*(k+1),2*k+1) = massmatrix(2*(k+1),2*k+1) + integral(@(x) f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% phi_left*psi_left
	massmatrix(2*(k+1),2*(k+1)+1) = massmatrix(2*(k+1),2*(k+1)+1) + integral(@(x) f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
    %psi_right*phi_right
    massmatrix(2*k+1,2*k) = massmatrix(2*k+1,2*k) + integral(@(x) f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));
    
    %psi_right*phi_left
    massmatrix(2*k+1,2*(k+1))=massmatrix(2*k+1,2*(k+1)) + integral(@(x) f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));
    
    
    
    %psi_left*phi_right
    massmatrix(2*(k+1)+1,2*k)=massmatrix(2*(k+1)+1,2*k) + integral(@(x) f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));
    
    %psi_left*psi_left
    massmatrix(2*(k+1)+1,2*(k+1))=massmatrix(2*(k+1)+1,2*(k+1)) + integral(@(x) f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));  
    
    
    
	%% psi_right*psi_right
	massmatrix(2*k+1,2*k+1) = massmatrix(2*k+1,2*k+1) + integral(@(x) f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% psi_right*psi_left
	massmatrix(2*k+1,2*(k+1)+1) = massmatrix(2*k+1,2*(k+1)+1) + integral(@(x) f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% psi_left*psi_right
	massmatrix(2*(k+1)+1,2*k+1) = massmatrix(2*(k+1)+1,2*k+1) + integral(@(x) f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% psi_left*psi_right
	massmatrix(2*(k+1)+1,2*(k+1)+1) = massmatrix(2*(k+1)+1,2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% Phi(k-1) n
	massrhs(2*k) =      massrhs(2*k) +        integral(@(x) f1((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));
	%% Phi(k) n
	massrhs(2*(k+1)) =  massrhs(2*(k+1)) +    integral(@(x) f2((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));
	%% Psi(k-1) n
	massrhs(2*k+1) =    massrhs(2*k+1) +      integral(@(x) f3((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));
	%% Psi(k) n
	massrhs(2*(k+1)+1) = massrhs(2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));

end

%Handle 1st interval 'half case'
massmatrix(1,1) = integral(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*(f1((x-node(1))) + alpha*f3((x-node(1)))),node(1),node(1)+dx);
massmatrix(1,2) = integral(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx);
massmatrix(2,1) = integral(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx);
massmatrix(1,3) = integral(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx);
massmatrix(3,1) = integral(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx);

massmatrix(2,2) = massmatrix(2,2) + integral(@(x) f2(x-node(1)).*f2(x-node(1)),node(1),node(1)+dx);
massmatrix(2,3) = massmatrix(2,3) + integral(@(x) f2(x-node(1)).*f4(x-node(1)),node(1),node(1)+dx);
massmatrix(3,2) = massmatrix(3,2) + integral(@(x) f2(x-node(1)).*f4(x-node(1)),node(1),node(1)+dx);
massmatrix(3,3) = massmatrix(3,3) + integral(@(x) f4(x-node(1)).*f4(x-node(1)),node(1),node(1)+dx);

massrhs(1) = massrhs(1) + integral(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*nGuess(x),node(1),node(1)+dx);
massrhs(2) = massrhs(2) + integral(@(x) f2((x-node(1))).*nGuess(x),node(1),node(1)+dx);
massrhs(3) = massrhs(3) + integral(@(x) f4((x-node(1))).*nGuess(x),node(1),node(1)+dx);

%Handle last interval 'half case'
last = 2*(numIntervals+1)-2;

massmatrix(last,last) = integral(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))),node(end)-dx,node(end));
massmatrix(last,last-2) = integral(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-2,last) = integral(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last,last-1) = integral(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-1,last) = integral(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end));

massmatrix(last-2,last-2) = massmatrix(last-2,last-2) + integral(@(x) f1((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-2,last-1) = massmatrix(last-2,last-1) + integral(@(x) f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-1,last-2) = massmatrix(last-1,last-2) + integral(@(x) f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-1,last-1) = massmatrix(last-1,last-1) + integral(@(x) f3((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end));

massrhs(last) = massrhs(last)+ integral(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*nGuess(x),node(end)-dx,node(end));

massrhs(last-2) = massrhs(last-2) + integral(@(x) f1((x-node(end)+dx)).*nGuess(x),node(end)-dx,node(end));
massrhs(last-1) = massrhs(last-1) + integral(@(x) f3((x-node(end)+dx)).*nGuess(x),node(end)-dx,node(end));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate composition of p(x) in finite element space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%alpha is boundary conditions term, p'(0)=alphap*p(0)
alphap = pprime(0)/pGuess(0);

%beta is boundary condition term, p'(L) = betap p(L)
betap= pprime(totalLength)/pGuess(totalLength);

massmatrixp = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

massrhsp = zeros(2*(numIntervals+1)-2,1);

for k = 1:numIntervals-2
	%%phi_righ*phi_right
	massmatrixp(2*k,2*k) = massmatrixp(2*k,2*k) + integral(@(x) f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_right*phi_left
	massmatrixp(2*k,2*(k+1)) = massmatrixp(2*k,2*(k+1)) + integral(@(x) f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

       
	%% phi_left*phi_right
	massmatrixp(2*(k+1),2*(k)) = massmatrixp(2*(k+1),2*(k)) + integral(@(x) f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*phi_left
	massmatrixp(2*(k+1),2*(k+1)) = massmatrixp(2*(k+1),2*(k+1)) + integral(@(x) f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError); 

    
    
	%% phi_right*psi_right
	massmatrixp(2*k,2*k+1) = massmatrixp(2*k,2*k+1) + integral(@(x) f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_right*psi_left
	massmatrixp(2*k,2*(k+1)+1) = massmatrixp(2*k,2*(k+1)+1) + integral(@(x) f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% phi_left*psi_right
	massmatrixp(2*(k+1),2*k+1) = massmatrixp(2*(k+1),2*k+1) + integral(@(x) f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*psi_left
	massmatrixp(2*(k+1),2*(k+1)+1) = massmatrixp(2*(k+1),2*(k+1)+1) + integral(@(x) f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
    %psi_right*phi_right
    massmatrixp(2*k+1,2*k) = massmatrixp(2*k+1,2*k) + integral(@(x) f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_right*phi_left
    massmatrixp(2*k+1,2*(k+1))=massmatrixp(2*k+1,2*(k+1)) + integral(@(x) f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    
    %psi_left*phi_right
    massmatrixp(2*(k+1)+1,2*k)=massmatrixp(2*(k+1)+1,2*k) + integral(@(x) f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_left*psi_left
    massmatrixp(2*(k+1)+1,2*(k+1))=massmatrixp(2*(k+1)+1,2*(k+1)) + integral(@(x) f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);  
    
    
    
	%% psi_right*psi_right
	massmatrixp(2*k+1,2*k+1) = massmatrixp(2*k+1,2*k+1) + integral(@(x) f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_right*psi_left
	massmatrixp(2*k+1,2*(k+1)+1) = massmatrixp(2*k+1,2*(k+1)+1) + integral(@(x) f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% psi_left*psi_right
	massmatrixp(2*(k+1)+1,2*k+1) = massmatrixp(2*(k+1)+1,2*k+1) + integral(@(x) f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_left*psi_right
	massmatrixp(2*(k+1)+1,2*(k+1)+1) = massmatrixp(2*(k+1)+1,2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% Phi(k-1) n
	massrhsp(2*k) =      massrhsp(2*k) +        integral(@(x) f1((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Phi(k) n
	massrhsp(2*(k+1)) =  massrhsp(2*(k+1)) +    integral(@(x) f2((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Psi(k-1) n
	massrhsp(2*k+1) =    massrhsp(2*k+1) +      integral(@(x) f3((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Psi(k) n
	massrhsp(2*(k+1)+1) = massrhsp(2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*pGuess(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

end
    
%Handle 1st interval 'half case'
massmatrixp(1,1) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(1,2) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(2,1) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(1,3) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(3,1) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

massmatrixp(2,2) = massmatrixp(2,2) + integral(@(x) f2(x).*f2(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(2,3) = massmatrixp(2,3) + integral(@(x) f2(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(3,2) = massmatrixp(3,2) + integral(@(x) f2(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(3,3) = massmatrixp(3,3) + integral(@(x) f4(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

massrhsp(1) = massrhsp(1) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*pGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massrhsp(2) = massrhsp(2) + integral(@(x) f2((x-node(1))).*pGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
massrhsp(3) = massrhsp(3) + integral(@(x) f4((x-node(1))).*pGuess(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
last = 2*(numIntervals+1)-2;
massmatrixp(last,last)   = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(last,last-2) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(last-2,last) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(last,last-1) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(last-1,last) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

massmatrixp(last-2,last-2) = massmatrixp(last-2,last-2) + integral(@(x) f1((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(last-2,last-1) = massmatrixp(last-2,last-1) + integral(@(x) f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(last-1,last-2) = massmatrixp(last-1,last-2) + integral(@(x) f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massmatrixp(last-1,last-1) = massmatrixp(last-1,last-1) + integral(@(x) f3((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

massrhsp(last) =   massrhsp(last)+    integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*pGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massrhsp(last-2) = massrhsp(last-2) + integral(@(x) f1((x-node(end)+dx)).*pGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
massrhsp(last-1) = massrhsp(last-1) +  integral(@(x) f3((x-node(end)+dx)).*pGuess(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform final calcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill scalling matrix
D = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

for index = 2:2:length(D)-1
   D(index,index) = 1;
   D(index+1,index+1) = 10^4*numIntervals;
end
D(1,1)=1;
D(length(D),length(D))=1;

%Scale Matrices using scaling matrix
LHSnscaled = D*massmatrix*D;
RHSnscaled = D*massrhs;

LHSpscaled = D*massmatrixp*D;
RHSpscaled = D*massrhsp;

ScaledCoeffsn = LHSnscaled\RHSnscaled;
Coeffsn = D'*ScaledCoeffsn;

ScaledCoeffsp = LHSpscaled\RHSpscaled;
Coeffsp = D'*ScaledCoeffsp;

%Solve simultaneous equations for basis function coefficients
%Coeffsn = massmatrix\massrhs;
%Coeffsp = massmatrixp\massrhsp;

%Collect Cj's in an array for display purposes
Cnjs = zeros(1,numIntervals+1);
Dnjs = zeros(1,numIntervals+1);

for index=2:length(Cnjs)-1
   Cnjs(1,index) = Coeffsn(2*index-2);
   Dnjs(1,index) = Coeffsn(2*index-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add in BC values omitted from calcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cnjs(1,1) = Coeffsn(1);
Dnjs(1,1) = alpha*Coeffsn(1);
Cnjs(1,numIntervals+1) = Coeffsn(length(Coeffsn));
Dnjs(1,numIntervals+1) = beta*Coeffsn(length(Coeffsn));


%Collect Cj's in an array for display purposes
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);

for index=2:length(Cnjs)-1
   Cpjs(1,index) = Coeffsp(2*index-2);
   Dpjs(1,index) = Coeffsp(2*index-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add in BC values omitted from calcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cpjs(1,1) = Coeffsp(1);
Dpjs(1,1) = alphap*Coeffsp(1);
Cpjs(1,numIntervals+1) = Coeffsp(length(Coeffsp));
Dpjs(1,numIntervals+1) = betap*Coeffsp(length(Coeffsp));

InitialCnjs = Cnjs;
InitialDnjs = Dnjs;
InitialCpjs = Cpjs;
InitialDpjs = Dpjs;

% fig4 = figure;
% plot(indexn,Cnjs,node,nGuess(node)),title('nGuess (recomposition)');
% saveas(fig4,strcat(filepath,'/nGuess_recomposition.jpg'))
% 
% fig4 = figure;
% plot(indexn,Dnjs,node,nprime(node)),title('nprime (recomposition)');
% saveas(fig4,strcat(filepath,'/nprime_recomposition.jpg'))
% 
% fig5 = figure;
% plot(indexn,Cpjs,node,pGuess(node)),title('pGuess (recomposition)');
% saveas(fig5,strcat(filepath,'/pGuess_recomposition.jpg'))
% 
% fig5 = figure;
% plot(indexn,Dpjs,node,pprime(node)),title('pprime (recomposition)');
% saveas(fig5,strcat(filepath,'/pprime_recomposition.jpg'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General BCs for calculated solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%alpha is boundary conditions term, n'(0)=alpha*n(0)
alphan = nprime(0)/nGuess(0);%Sr/Dn;

%beta is boundary condition term, n'(L) = beta n(L)
betan = nprime(totalLength)/nGuess(totalLength);%-Sr/Dn;

%alpha is boundary conditions term, p'(0)=alpha*p(0)
alphap = pprime(0)/pGuess(0);%Sr/Dp;
%beta is boundary condition term, p'(L) = beta p(L)
betap = pprime(totalLength)/pGuess(totalLength);%-Sr/Dp;

%pGuess = matlabFunction(pGuess_syms);
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

for k = 1:numIntervals-2
   	%% Phi(k-1) n
	rhs_matrix(2*k) =      rhs_matrix(2*k) +        integral(@(x) f1((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Phi(k) n
	rhs_matrix(2*(k+1)) =  rhs_matrix(2*(k+1)) +    integral(@(x) f2((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Psi(k-1) n
	rhs_matrix(2*k+1) =    rhs_matrix(2*k+1) +      integral(@(x) f3((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Psi(k) n
	rhs_matrix(2*(k+1)+1) = rhs_matrix(2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);    
end

rhs_matrix(1) = rhs_matrix(1) + integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*RHS(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(2) = rhs_matrix(2) + integral(@(x) f2((x-node(1))).*RHS(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(3) = rhs_matrix(3) + integral(@(x) f4((x-node(1))).*RHS(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

rhs_matrix(last-2) = rhs_matrix(last-2) + integral(@(x) f1((x-node(end)+dx)).*RHS(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last-1) = rhs_matrix(last-1) + integral(@(x) f3((x-node(end)+dx)).*RHS(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last) = rhs_matrix(last) + integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*RHS(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
 
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

matrix1 = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);


for k = 1:numIntervals-2
	%%phi_righ*phi_right
	matrix1(2*k,2*k) = matrix1(2*k,2*k) + integral(@(x) f1((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_right*phi_left
	matrix1(2*k,2*(k+1)) = matrix1(2*k,2*(k+1)) + integral(@(x) f1((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

       
	%% phi_left*phi_right
	matrix1(2*(k+1),2*(k)) = matrix1(2*(k+1),2*(k)) + integral(@(x) f2((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*phi_left
	matrix1(2*(k+1),2*(k+1)) = matrix1(2*(k+1),2*(k+1)) + integral(@(x) f2((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError); 

    
    
	%% phi_right*psi_right
	matrix1(2*k,2*k+1) = matrix1(2*k,2*k+1) + integral(@(x) f1((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

	%% phi_right*psi_left
	matrix1(2*k,2*(k+1)+1) = matrix1(2*k,2*(k+1)+1) + integral(@(x) f1((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% phi_left*psi_right
	matrix1(2*(k+1),2*k+1) = matrix1(2*(k+1),2*k+1) + integral(@(x) f2((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*psi_left
	matrix1(2*(k+1),2*(k+1)+1) = matrix1(2*(k+1),2*(k+1)+1) + integral(@(x) f2((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

    
    
    %psi_right*phi_right
    matrix1(2*k+1,2*k) = matrix1(2*k+1,2*k) + integral(@(x)  f3((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%
    
    %psi_right*phi_left
    matrix1(2*k+1,2*(k+1))=matrix1(2*k+1,2*(k+1)) + integral(@(x) f3((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    
    %psi_left*phi_right
    matrix1(2*(k+1)+1,2*k)=matrix1(2*(k+1)+1,2*k) + integral(@(x) f4((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_left*psi_left
    matrix1(2*(k+1)+1,2*(k+1))=matrix1(2*(k+1)+1,2*(k+1)) + integral(@(x) f4((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);  %
    
    
    
	%% psi_right*psi_right
	matrix1(2*k+1,2*k+1) = matrix1(2*k+1,2*k+1) + integral(@(x) f3((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_right*psi_left
	matrix1(2*k+1,2*(k+1)+1) = matrix1(2*k+1,2*(k+1)+1) + integral(@(x) f3((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% psi_left*psi_right
	matrix1(2*(k+1)+1,2*k+1) = matrix1(2*(k+1)+1,2*k+1) + integral(@(x) f4((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_left*psi_right
	matrix1(2*(k+1)+1,2*(k+1)+1) = matrix1(2*(k+1)+1,2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end

%Handle 1st interval 'half case'
matrix1(1,1) = integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*(f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(1,2) = integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*f2primeprime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(1,3) = integral(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*f4primeprime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(2,1) = integral(@(x) (f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(3,1) = integral(@(x) (f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrix1(2,2) = matrix1(2,2) + integral(@(x) f2(x).*f2primeprime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(2,3) = matrix1(2,3) + integral(@(x) f2(x).*f4primeprime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(3,2) = matrix1(3,2) + integral(@(x) f2primeprime(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(3,3) = matrix1(3,3) + integral(@(x) f4(x).*f4primeprime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
matrix1(last,last) =   integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*(f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last,last-2) = integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1primeprime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last,last-1) = integral(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3primeprime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-2,last) = integral(@(x) (f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-1,last) = integral(@(x) (f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrix1(last-2,last-2) = matrix1(last-2,last-2) + integral(@(x) f1primeprime((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-2,last-1) = matrix1(last-2,last-1) + integral(@(x) f1((x-node(end)+dx)).*f3primeprime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-1,last-2) = matrix1(last-1,last-2) + integral(@(x) f1primeprime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-1,last-1) = matrix1(last-1,last-1) + integral(@(x) f3primeprime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);   

   Matrix1Finaln = Cn1*matrix1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End matrix of n''(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates for RHS of equation (forcing function times test functions)
%This calculation is performed outside the iterative loop, as it will not
%change over the course of changing p(x), n(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RHS = Gopp;

%Generate matrix of test functions against RHS of equation
rhs_matrix = zeros(2*(numIntervals+1)-2,1);

for k = 1:numIntervals-2
   	%% Phi(k-1) n
	rhs_matrix(2*k) =      rhs_matrix(2*k) +        integral(@(x) f1((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Phi(k) n
	rhs_matrix(2*(k+1)) =  rhs_matrix(2*(k+1)) +    integral(@(x) f2((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Psi(k-1) n
	rhs_matrix(2*k+1) =    rhs_matrix(2*k+1) +      integral(@(x) f3((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
	%% Psi(k) n
	rhs_matrix(2*(k+1)+1) = rhs_matrix(2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*RHS(x),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);    
end

rhs_matrix(1) = rhs_matrix(1) + integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*RHS(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(2) = rhs_matrix(2) + integral(@(x) f2((x-node(1))).*RHS(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(3) = rhs_matrix(3) + integral(@(x) f4((x-node(1))).*RHS(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

rhs_matrix(last-2) = rhs_matrix(last-2) + integral(@(x) f1((x-node(end)+dx)).*RHS(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last-1) = rhs_matrix(last-1) + integral(@(x) f3((x-node(end)+dx)).*RHS(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
rhs_matrix(last) = rhs_matrix(last) + integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*RHS(x),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
 
RHSOutp =rhs_matrix;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End RHS calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Fill matrix associated with p''(x) term
% %Performed outside iterative loop because this term is not influenced by
% %p(x) or n(x).  Should not change over calculations.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Front end coefficient of n''(x) term
    Cp1=ksi*T*mup/q-Dp;

matrix1 = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);


for k = 1:numIntervals-2
	%%phi_righ*phi_right
	matrix1(2*k,2*k) = matrix1(2*k,2*k) + integral(@(x) f1((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_right*phi_left
	matrix1(2*k,2*(k+1)) = matrix1(2*k,2*(k+1)) + integral(@(x) f1((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

       
	%% phi_left*phi_right
	matrix1(2*(k+1),2*(k)) = matrix1(2*(k+1),2*(k)) + integral(@(x) f2((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*phi_left
	matrix1(2*(k+1),2*(k+1)) = matrix1(2*(k+1),2*(k+1)) + integral(@(x) f2((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError); 

    
    
	%% phi_right*psi_right
	matrix1(2*k,2*k+1) = matrix1(2*k,2*k+1) + integral(@(x) f1((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

	%% phi_right*psi_left
	matrix1(2*k,2*(k+1)+1) = matrix1(2*k,2*(k+1)+1) + integral(@(x) f1((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% phi_left*psi_right
	matrix1(2*(k+1),2*k+1) = matrix1(2*(k+1),2*k+1) + integral(@(x) f2((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*psi_left
	matrix1(2*(k+1),2*(k+1)+1) = matrix1(2*(k+1),2*(k+1)+1) + integral(@(x) f2((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

    
    
    %psi_right*phi_right
    matrix1(2*k+1,2*k) = matrix1(2*k+1,2*k) + integral(@(x)  f3((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%
    
    %psi_right*phi_left
    matrix1(2*k+1,2*(k+1))=matrix1(2*k+1,2*(k+1)) + integral(@(x) f3((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    
    %psi_left*phi_right
    matrix1(2*(k+1)+1,2*k)=matrix1(2*(k+1)+1,2*k) + integral(@(x) f4((x-node(k+1))).*f1primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_left*psi_left
    matrix1(2*(k+1)+1,2*(k+1))=matrix1(2*(k+1)+1,2*(k+1)) + integral(@(x) f4((x-node(k+1))).*f2primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);  %
    
    
    
	%% psi_right*psi_right
	matrix1(2*k+1,2*k+1) = matrix1(2*k+1,2*k+1) + integral(@(x) f3((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_right*psi_left
	matrix1(2*k+1,2*(k+1)+1) = matrix1(2*k+1,2*(k+1)+1) + integral(@(x) f3((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% psi_left*psi_right
	matrix1(2*(k+1)+1,2*k+1) = matrix1(2*(k+1)+1,2*k+1) + integral(@(x) f4((x-node(k+1))).*f3primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_left*psi_right
	matrix1(2*(k+1)+1,2*(k+1)+1) = matrix1(2*(k+1)+1,2*(k+1)+1) + integral(@(x) f4((x-node(k+1))).*f4primeprime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end

%Handle 1st interval 'half case'
matrix1(1,1) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*(f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(1,2) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f2primeprime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(1,3) = integral(@(x) (f1((x-node(1))) + alphap*f3((x-node(1)))).*f4primeprime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(2,1) = integral(@(x) (f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(3,1) = integral(@(x) (f1primeprime((x-node(1))) + alphap*f3primeprime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrix1(2,2) = matrix1(2,2) + integral(@(x) f2(x).*f2primeprime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(2,3) = matrix1(2,3) + integral(@(x) f2(x).*f4primeprime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(3,2) = matrix1(3,2) + integral(@(x) f2primeprime(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix1(3,3) = matrix1(3,3) + integral(@(x) f4(x).*f4primeprime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

%Handle last interval 'half case'
matrix1(last,last) =   integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*(f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last,last-2) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1primeprime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last,last-1) = integral(@(x) (f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3primeprime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-2,last) = integral(@(x) (f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-1,last) = integral(@(x) (f2primeprime((x-node(end)+dx)) + betap*f4primeprime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrix1(last-2,last-2) = matrix1(last-2,last-2) + integral(@(x) f1primeprime((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-2,last-1) = matrix1(last-2,last-1) + integral(@(x) f1((x-node(end)+dx)).*f3primeprime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-1,last-2) = matrix1(last-1,last-2) + integral(@(x) f1primeprime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix1(last-1,last-1) = matrix1(last-1,last-1) + integral(@(x) f3primeprime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);   

    Matrix1Finalp = Cp1*matrix1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of p''(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vestiges from more expansive versions of code
%Output matrices from these sections should empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterationCount=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate current representations of p(x),p'(x),p''(x)
%Values associated with 1st node in system are undefined
%Value defaults to 0, but is never used in code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_dist = cell(1,numIntervals);
% pprime_dist =cell(1,numIntervals);
% pdouble_dist =cell(1,numIntervals);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Finds functions in the interval to the left of the kth node
% %Note: This implies that p(1), p'(1), p''(1) are undefined
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for index=1:numIntervals
% p_dist{index}= px_v2(index+1,node,Cpjs,Dpjs,f1,f2,f3,f4);
% pprime_dist{index}=px_v2(index+1,node,Cpjs,Dpjs,f1prime,f2prime,f3prime,f4prime);
% pdouble_dist{index}=px_v2(index+1,node,Cpjs,Dpjs,f1primeprime,f2primeprime,f3primeprime,f4primeprime);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with n'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iterationCount=1:maxIterations
%Front end coefficient for this matrix.  DOES NOT INCLUDE CONTRIBUTION FROM
%p(x) DISTRIBUTION.
     Cn2 = ksi*T*mun/q;

   matrix2 = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%pprime_dist{k}(x)./p_dist{k}(x).*
for k = 1:numIntervals-2
    k=k+1;
    %Define p'(x)/p(x) for this interval
    pprime_over_p = @(x) (Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
	k=k-1;
    
    %%phi_righ*phi_right
	matrix2(2*k,2*k) = matrix2(2*k,2*k) + integral(@(x) pprime_over_p(x).*f1((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
	%% phi_right*phi_left
	matrix2(2*k,2*(k+1)) = matrix2(2*k,2*(k+1)) + integral(@(x) pprime_over_p(x).*f1((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

       
	%% phi_left*phi_right
	matrix2(2*(k+1),2*(k)) = matrix2(2*(k+1),2*(k)) + integral(@(x) pprime_over_p(x).*f2((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*phi_left
	matrix2(2*(k+1),2*(k+1)) = matrix2(2*(k+1),2*(k+1)) + integral(@(x) pprime_over_p(x).*f2((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError); 

    
    
	%% phi_right*psi_right
	matrix2(2*k,2*k+1) = matrix2(2*k,2*k+1) + integral(@(x) pprime_over_p(x).*f1((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

	%% phi_right*psi_left
	matrix2(2*k,2*(k+1)+1) = matrix2(2*k,2*(k+1)+1) + integral(@(x) pprime_over_p(x).*f1((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% phi_left*psi_right
	matrix2(2*(k+1),2*k+1) = matrix2(2*(k+1),2*k+1) + integral(@(x) pprime_over_p(x).*f2((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*psi_left
	matrix2(2*(k+1),2*(k+1)+1) = matrix2(2*(k+1),2*(k+1)+1) + integral(@(x) pprime_over_p(x).*f2((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

    
    
    %psi_right*phi_right
    matrix2(2*k+1,2*k) = matrix2(2*k+1,2*k) + integral(@(x)  pprime_over_p(x).*f3((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%
    
    %psi_right*phi_left
    matrix2(2*k+1,2*(k+1))=matrix2(2*k+1,2*(k+1)) + integral(@(x) pprime_over_p(x).*f3((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    
    %psi_left*phi_right
    matrix2(2*(k+1)+1,2*k)=matrix2(2*(k+1)+1,2*k) + integral(@(x) pprime_over_p(x).*f4((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_left*psi_left
    matrix2(2*(k+1)+1,2*(k+1))=matrix2(2*(k+1)+1,2*(k+1)) + integral(@(x) pprime_over_p(x).*f4((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);  %
    
    
    
	%% psi_right*psi_right
	matrix2(2*k+1,2*k+1) = matrix2(2*k+1,2*k+1) + integral(@(x) pprime_over_p(x).*f3((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_right*psi_left
	matrix2(2*k+1,2*(k+1)+1) = matrix2(2*k+1,2*(k+1)+1) + integral(@(x) pprime_over_p(x).*f3((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% psi_left*psi_right
	matrix2(2*(k+1)+1,2*k+1) = matrix2(2*(k+1)+1,2*k+1) + integral(@(x) pprime_over_p(x).*f4((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_left*psi_right
	matrix2(2*(k+1)+1,2*(k+1)+1) = matrix2(2*(k+1)+1,2*(k+1)+1) + integral(@(x) pprime_over_p(x).*f4((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end
k=1;
%Define p'(x)/p(x) for this interval
pprime_over_p = @(x) (Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));
    
%Handle 1st interval 'half case'
matrix2(1,1) = integral(@(x) pprime_over_p(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(1,2) = integral(@(x) pprime_over_p(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f2prime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(1,3) = integral(@(x) pprime_over_p(x).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f4prime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(2,1) = integral(@(x) pprime_over_p(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(3,1) = integral(@(x) pprime_over_p(x).*(f1prime((x-node(1))) + alphan*f3prime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrix2(2,2) = matrix2(2,2) + integral(@(x) pprime_over_p(x).*f2(x).*f2prime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(2,3) = matrix2(2,3) + integral(@(x) pprime_over_p(x).*f2(x).*f4prime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(3,2) = matrix2(3,2) + integral(@(x) pprime_over_p(x).*f2prime(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(3,3) = matrix2(3,3) + integral(@(x) pprime_over_p(x).*f4(x).*f4prime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

k=numIntervals;
%Define p'(x)/p(x) for this interval
pprime_over_p = @(x) (Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)));	
    
%Handle last interval 'half case'
matrix2(last,last) =   integral(@(x) pprime_over_p(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last,last-2) = integral(@(x) pprime_over_p(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1prime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last,last-1) = integral(@(x) pprime_over_p(x).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3prime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-2,last) = integral(@(x) pprime_over_p(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-1,last) = integral(@(x) pprime_over_p(x).*(f2prime((x-node(end)+dx)) + betan*f4prime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrix2(last-2,last-2) = matrix2(last-2,last-2) + integral(@(x) pprime_over_p(x).*f1prime((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-2,last-1) = matrix2(last-2,last-1) + integral(@(x) pprime_over_p(x).*f1((x-node(end)+dx)).*f3prime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-1,last-2) = matrix2(last-1,last-2) + integral(@(x) pprime_over_p(x).*f1prime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-1,last-1) = matrix2(last-1,last-1) + integral(@(x) pprime_over_p(x).* f3prime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);   


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
 matrix3 = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

for k = 1:numIntervals-2
    k=k+1;
    %Define p'(x)/p(x) for this interval
    pBlock = @(x) ((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))-((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2);
	k=k-1;
    
	%%phi_righ*phi_right
	matrix3(2*k,2*k) = matrix3(2*k,2*k) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
	%% phi_right*phi_left
	matrix3(2*k,2*(k+1)) = matrix3(2*k,2*(k+1)) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

       
	%% phi_left*phi_right
	matrix3(2*(k+1),2*(k)) = matrix3(2*(k+1),2*(k)) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*phi_left
	matrix3(2*(k+1),2*(k+1)) = matrix3(2*(k+1),2*(k+1)) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError); 

    
    
	%% phi_right*psi_right
	matrix3(2*k,2*k+1) = matrix3(2*k,2*k+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

	%% phi_right*psi_left
	matrix3(2*k,2*(k+1)+1) = matrix3(2*k,2*(k+1)+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% phi_left*psi_right
	matrix3(2*(k+1),2*k+1) = matrix3(2*(k+1),2*k+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*psi_left
	matrix3(2*(k+1),2*(k+1)+1) = matrix3(2*(k+1),2*(k+1)+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

    
    
    %psi_right*phi_right
    matrix3(2*k+1,2*k) = matrix3(2*k+1,2*k) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%
    
    %psi_right*phi_left
    matrix3(2*k+1,2*(k+1))=matrix3(2*k+1,2*(k+1)) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    
    %psi_left*phi_right
    matrix3(2*(k+1)+1,2*k)=matrix3(2*(k+1)+1,2*k) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_left*psi_left
    matrix3(2*(k+1)+1,2*(k+1))=matrix3(2*(k+1)+1,2*(k+1)) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);  %
    
    
    
	%% psi_right*psi_right
	matrix3(2*k+1,2*k+1) = matrix3(2*k+1,2*k+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_right*psi_left
	matrix3(2*k+1,2*(k+1)+1) = matrix3(2*k+1,2*(k+1)+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% psi_left*psi_right
	matrix3(2*(k+1)+1,2*k+1) = matrix3(2*(k+1)+1,2*k+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_left*psi_right
	matrix3(2*(k+1)+1,2*(k+1)+1) = matrix3(2*(k+1)+1,2*(k+1)+1) + integral(@(x) (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end

k=1;
%Define p'(x)/p(x) for this interval
pBlock = @(x) ((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))-((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2);
	

%Handle 1st interval 'half case'
matrix3(1,1) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(1,2) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(1,3) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(2,1) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(3,1) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f1((x-node(1))) + alphan*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

 matrix3(2,2) = matrix3(2,2) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f2(x).*f2(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
 matrix3(2,3) = matrix3(2,3) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f2(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
 matrix3(3,2) = matrix3(3,2) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f2(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
 matrix3(3,3) = matrix3(3,3) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f4(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

k=numIntervals;
%Define p'(x)/p(x) for this interval
pBlock = @(x) ((Cpjs(k).*f1primeprime(x-node(k)) + Cpjs(k+1).*f2primeprime(x-node(k)) + Dpjs(k).*f3primeprime(x-node(k)) + Dpjs(k+1).*f4primeprime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))-((Cpjs(k).*f1prime(x-node(k)) + Cpjs(k+1).*f2prime(x-node(k)) + Dpjs(k).*f3prime(x-node(k)) + Dpjs(k+1).*f4prime(x-node(k)))./(Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k)))).^2);
	
%Handle last interval 'half case'
matrix3(last,last) =   integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last,last-2) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last,last-1) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-2,last) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-1,last) = integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrix3(last-2,last-2) = matrix3(last-2,last-2) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f1((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-2,last-1) = matrix3(last-2,last-1) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-1,last-2) = matrix3(last-1,last-2) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-1,last-1) = matrix3(last-1,last-1) + integral(@(x)  (ksi.*T.*mun./q.*(pBlock(x))-1/taun).*f3((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);   


Matrix3Finaln = matrix3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of n(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%
% 
%Combine expressions from each portion for total matrix of LHS
LHSn = Matrix1Finaln +Matrix2Finaln +Matrix3Finaln;

%Scale Matrices using scaling matrix
LHSscaled = D*LHSn*D;
RHSscaled = D*RHSOutn;

ScaledCoeffs = LHSscaled\RHSscaled;
Coeffsn = D'*ScaledCoeffs;
%Solve simultaneous equations for basis function coefficients
%Coeffsn = LHSn\RHSOutn;

%Collect Cj's in an array for display purposes
Cnjs = zeros(1,numIntervals+1);
Dnjs = zeros(1,numIntervals+1);

for index=2:length(Cnjs)-1
   Cnjs(1,index) = Coeffsn(2*index-2);
   Dnjs(1,index) = Coeffsn(2*index-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add in BC values omitted from calcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cnjs(1,1) = Coeffsn(1);
Dnjs(1,1) = alphan*Coeffsn(1);
Cnjs(1,numIntervals+1) = Coeffsn(length(Coeffsn));
Dnjs(1,numIntervals+1) = betan*Coeffsn(length(Coeffsn));


fig4 = figure;
plot(indexn,Cnjs,'x',indexn,nGuess(indexn));title('n(x) calculated'),legend('n(x) - calculated','n(x) - original');set(gca,'FontSize',16),xlabel('Depth (in cm)'),ylabel('Carrier concentration (in cm^-3)');
saveas(fig4,strcat(filepath,'/n_calculated_iteration',num2str(iterationCount),'.jpg'))

fig6 = figure;
plot(indexn,Dnjs,'x',indexn,nprime(indexn)),title('nprime calculated'),legend('Dnjs','nprime - original');
saveas(fig6,strcat(filepath,'/nprime_calculated_iteration',num2str(iterationCount),'.jpg'))

fig5 = figure;
plot_together(Cnjs,Dnjs,totalLength,numIntervals),title('plot together');
saveas(fig5,strcat(filepath,'/plot_together_n',num2str(iterationCount),'.jpg'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vestiges from more expansive versions of code
%Output matrices from these sections should empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate current representations of p(x),p'(x),p''(x)
%Values associated with 1st node in system are undefined
%Value defaults to 0, but is never used in code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_dist = cell(1,numIntervals);
% pprime_dist =cell(1,numIntervals);
% pdouble_dist =cell(1,numIntervals);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Finds functions in the interval to the left of the kth node
% %Note: This implies that p(1), p'(1), p''(1) are undefined
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for index=1:numIntervals
% p_dist{index}= px_v2(index+1,node,Cpjs,Dpjs,f1,f2,f3,f4);
% pprime_dist{index}=px_v2(index+1,node,Cpjs,Dpjs,f1prime,f2prime,f3prime,f4prime);
% pdouble_dist{index}=px_v2(index+1,node,Cpjs,Dpjs,f1primeprime,f2primeprime,f3primeprime,f4primeprime);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with n'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Front end coefficient for this matrix.  DOES NOT INCLUDE CONTRIBUTION FROM
%p(x) DISTRIBUTION.
     Cp2 = -ksi*T*mup/q;

   matrix2 = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

%pprime_dist{k}(x)./p_dist{k}(x).*
for k = 1:numIntervals-2
    k=k+1;
    %Define p'(x)/p(x) for this interval
    nprime_over_n = @(x) (Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
	k=k-1;
    
    %%phi_righ*phi_right
	matrix2(2*k,2*k) = matrix2(2*k,2*k) + integral(@(x) nprime_over_n(x).*f1((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
	%% phi_right*phi_left
	matrix2(2*k,2*(k+1)) = matrix2(2*k,2*(k+1)) + integral(@(x) nprime_over_n(x).*f1((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

       
	%% phi_left*phi_right
	matrix2(2*(k+1),2*(k)) = matrix2(2*(k+1),2*(k)) + integral(@(x) nprime_over_n(x).*f2((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*phi_left
	matrix2(2*(k+1),2*(k+1)) = matrix2(2*(k+1),2*(k+1)) + integral(@(x) nprime_over_n(x).*f2((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError); 

    
    
	%% phi_right*psi_right
	matrix2(2*k,2*k+1) = matrix2(2*k,2*k+1) + integral(@(x) nprime_over_n(x).*f1((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

	%% phi_right*psi_left
	matrix2(2*k,2*(k+1)+1) = matrix2(2*k,2*(k+1)+1) + integral(@(x) nprime_over_n(x).*f1((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% phi_left*psi_right
	matrix2(2*(k+1),2*k+1) = matrix2(2*(k+1),2*k+1) + integral(@(x) nprime_over_n(x).*f2((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*psi_left
	matrix2(2*(k+1),2*(k+1)+1) = matrix2(2*(k+1),2*(k+1)+1) + integral(@(x) nprime_over_n(x).*f2((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

    
    
    %psi_right*phi_right
    matrix2(2*k+1,2*k) = matrix2(2*k+1,2*k) + integral(@(x)  nprime_over_n(x).*f3((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%
    
    %psi_right*phi_left
    matrix2(2*k+1,2*(k+1))=matrix2(2*k+1,2*(k+1)) + integral(@(x) nprime_over_n(x).*f3((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    
    %psi_left*phi_right
    matrix2(2*(k+1)+1,2*k)=matrix2(2*(k+1)+1,2*k) + integral(@(x) nprime_over_n(x).*f4((x-node(k+1))).*f1prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_left*psi_left
    matrix2(2*(k+1)+1,2*(k+1))=matrix2(2*(k+1)+1,2*(k+1)) + integral(@(x) nprime_over_n(x).*f4((x-node(k+1))).*f2prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);  %
    
    
    
	%% psi_right*psi_right
	matrix2(2*k+1,2*k+1) = matrix2(2*k+1,2*k+1) + integral(@(x) nprime_over_n(x).*f3((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_right*psi_left
	matrix2(2*k+1,2*(k+1)+1) = matrix2(2*k+1,2*(k+1)+1) + integral(@(x) nprime_over_n(x).*f3((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% psi_left*psi_right
	matrix2(2*(k+1)+1,2*k+1) = matrix2(2*(k+1)+1,2*k+1) + integral(@(x) nprime_over_n(x).*f4((x-node(k+1))).*f3prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_left*psi_right
	matrix2(2*(k+1)+1,2*(k+1)+1) = matrix2(2*(k+1)+1,2*(k+1)+1) + integral(@(x) nprime_over_n(x).*f4((x-node(k+1))).*f4prime((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end
k=1;
%Define p'(x)/p(x) for this interval
nprime_over_n = @(x) (Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
	    
%Handle 1st interval 'half case'
matrix2(1,1) = integral(@(x) nprime_over_n(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(1,2) = integral(@(x) nprime_over_n(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f2prime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(1,3) = integral(@(x) nprime_over_n(x).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f4prime((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(2,1) = integral(@(x) nprime_over_n(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(3,1) = integral(@(x) nprime_over_n(x).*(f1prime((x-node(1))) + alphap*f3prime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

matrix2(2,2) = matrix2(2,2) + integral(@(x) nprime_over_n(x).*f2(x).*f2prime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(2,3) = matrix2(2,3) + integral(@(x) nprime_over_n(x).*f2(x).*f4prime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(3,2) = matrix2(3,2) + integral(@(x) nprime_over_n(x).*f2prime(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix2(3,3) = matrix2(3,3) + integral(@(x) nprime_over_n(x).*f4(x).*f4prime(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

k=numIntervals;
%Define p'(x)/p(x) for this interval
nprime_over_n = @(x) (Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)));
	   
%Handle last interval 'half case'
matrix2(last,last) =   integral(@(x) nprime_over_n(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last,last-2) = integral(@(x) nprime_over_n(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1prime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last,last-1) = integral(@(x) nprime_over_n(x).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3prime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-2,last) = integral(@(x) nprime_over_n(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-1,last) = integral(@(x) nprime_over_n(x).*(f2prime((x-node(end)+dx)) + betap*f4prime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrix2(last-2,last-2) = matrix2(last-2,last-2) + integral(@(x) nprime_over_n(x).*f1prime((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-2,last-1) = matrix2(last-2,last-1) + integral(@(x) nprime_over_n(x).*f1((x-node(end)+dx)).*f3prime((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-1,last-2) = matrix2(last-1,last-2) + integral(@(x) nprime_over_n(x).*f1prime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix2(last-1,last-1) = matrix2(last-1,last-1) + integral(@(x) nprime_over_n(x).* f3prime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);   


   Matrix2Finalp = Cp2*matrix2;
  
  
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
 matrix3 = sparse(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

for k = 1:numIntervals-2
    k=k+1;
    %Define p'(x)/p(x) for this interval
    nBlock = @(x) (((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
	k=k-1;
    
	%%phi_righ*phi_right
	matrix3(2*k,2*k) = matrix3(2*k,2*k) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
	%% phi_right*phi_left
	matrix3(2*k,2*(k+1)) = matrix3(2*k,2*(k+1)) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

       
	%% phi_left*phi_right
	matrix3(2*(k+1),2*(k)) = matrix3(2*(k+1),2*(k)) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*phi_left
	matrix3(2*(k+1),2*(k+1)) = matrix3(2*(k+1),2*(k+1)) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError); 

    
    
	%% phi_right*psi_right
	matrix3(2*k,2*k+1) = matrix3(2*k,2*k+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

	%% phi_right*psi_left
	matrix3(2*k,2*(k+1)+1) = matrix3(2*k,2*(k+1)+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% phi_left*psi_right
	matrix3(2*(k+1),2*k+1) = matrix3(2*(k+1),2*k+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% phi_left*psi_left
	matrix3(2*(k+1),2*(k+1)+1) = matrix3(2*(k+1),2*(k+1)+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%

    
    
    %psi_right*phi_right
    matrix3(2*k+1,2*k) = matrix3(2*k+1,2*k) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);%
    
    %psi_right*phi_left
    matrix3(2*k+1,2*(k+1))=matrix3(2*k+1,2*(k+1)) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    
    
    %psi_left*phi_right
    matrix3(2*(k+1)+1,2*k)=matrix3(2*(k+1)+1,2*k) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
    
    %psi_left*psi_left
    matrix3(2*(k+1)+1,2*(k+1))=matrix3(2*(k+1)+1,2*(k+1)) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);  %
    
    
    
	%% psi_right*psi_right
	matrix3(2*k+1,2*k+1) = matrix3(2*k+1,2*k+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_right*psi_left
	matrix3(2*k+1,2*(k+1)+1) = matrix3(2*k+1,2*(k+1)+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

    
    
	%% psi_left*psi_right
	matrix3(2*(k+1)+1,2*k+1) = matrix3(2*(k+1)+1,2*k+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);

	%% psi_left*psi_right
	matrix3(2*(k+1)+1,2*(k+1)+1) = matrix3(2*(k+1)+1,2*(k+1)+1) + integral(@(x) (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2),'AbsTol',AbsError,'RelTol',RelError);
end

k=1;
%Define p'(x)/p(x) for this interval
nBlock = @(x) (((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
		

%Handle 1st interval 'half case'
matrix3(1,1) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*(f1((x-node(1))) + alphap*f3((x-node(1)))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(1,2) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(1,3) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(2,1) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
matrix3(3,1) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f1((x-node(1))) + alphap*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

 matrix3(2,2) = matrix3(2,2) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f2(x).*f2(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
 matrix3(2,3) = matrix3(2,3) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f2(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
 matrix3(3,2) = matrix3(3,2) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f2(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);
 matrix3(3,3) = matrix3(3,3) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f4(x).*f4(x),node(1),node(1)+dx,'AbsTol',AbsError,'RelTol',RelError);

k=numIntervals;
%Define p'(x)/p(x) for this interval
nBlock = @(x) (((Cnjs(k).*f1prime(x-node(k)) + Cnjs(k+1).*f2prime(x-node(k)) + Dnjs(k).*f3prime(x-node(k)) + Dnjs(k+1).*f4prime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k)))).^2-(Cnjs(k).*f1primeprime(x-node(k)) + Cnjs(k+1).*f2primeprime(x-node(k)) + Dnjs(k).*f3primeprime(x-node(k)) + Dnjs(k+1).*f4primeprime(x-node(k)))./(Cnjs(k).*f1(x-node(k)) + Cnjs(k+1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k)) + Dnjs(k+1).*f4(x-node(k))));
		
%Handle last interval 'half case'
matrix3(last,last) =   integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last,last-2) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last,last-1) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-2,last) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-1,last) = integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*(f2((x-node(end)+dx)) + betap*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);

matrix3(last-2,last-2) = matrix3(last-2,last-2) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f1((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-2,last-1) = matrix3(last-2,last-1) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-1,last-2) = matrix3(last-1,last-2) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);
matrix3(last-1,last-1) = matrix3(last-1,last-1) + integral(@(x)  (ksi.*T.*mup./q.*(nBlock(x))-1/taun).*f3((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end),'AbsTol',AbsError,'RelTol',RelError);   


Matrix3Finalp = matrix3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of n(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%
% 
%Combine expressions from each portion for total matrix of LHS
LHSp = Matrix1Finalp +Matrix2Finalp +Matrix3Finalp;

%Scale Matrices using scaling matrix
LHSscaled = D*LHSp*D;
RHSscaled = D*RHSOutp;

ScaledCoeffs = LHSscaled\RHSscaled;
Coeffsp = D'*ScaledCoeffs;

%Solve simultaneous equations for basis function coefficients
%Coeffsp = LHSp\RHSOutp;

%Collect Cj's in an array for display purposes
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);

for index=2:length(Cnjs)-1
   Cpjs(1,index) = Coeffsp(2*index-2);
   Dpjs(1,index) = Coeffsp(2*index-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add in BC values omitted from calcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cpjs(1,1) = Coeffsp(1);
Dpjs(1,1) = alphap*Coeffsp(1);
Cpjs(1,numIntervals+1) = Coeffsp(length(Coeffsn));
Dpjs(1,numIntervals+1) = betap*Coeffsp(length(Coeffsn));
% 

fig4 = figure;
plot(indexn,Cpjs,'x',indexn,pGuess(indexn),indexn,pOrig(indexn));title('p(x) calculated'),legend('p(x) - calculated','p(x) - input','p(x) - original');set(gca,'FontSize',16),xlabel('Depth (in cm)'),ylabel('Carrier concentration (in cm^-3)');
saveas(fig4,strcat(filepath,'/p_calculated_iteration',num2str(iterationCount),'.jpg'))

fig6 = figure;
plot(indexn,Dpjs,'x',indexn,pprime(indexn),indexn,pprimeOrig(indexn)),title('pprime calculated'),legend('Dpjs','pGuess - error','pGuess - orignial');
saveas(fig6,strcat(filepath,'/pprime_calculated_iteration',num2str(iterationCount),'.jpg'))

fig5 = figure;
plot_together(Cnjs,Dnjs,totalLength,numIntervals),title('plot together');
saveas(fig5,strcat(filepath,'/plot_together_p',num2str(iterationCount),'.jpg'))

save(strcat(filepath,'/','iteration',num2str(iterationCount),'.mat'))
end

%Displays runtime of m-file
time = cputime-t



%Store variable values to be used later
save(strcat(filepath,'/',num2str(numIntervals),'terms.mat'))