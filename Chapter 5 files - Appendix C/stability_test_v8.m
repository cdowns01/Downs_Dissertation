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
dx = intervalLength;
maxIterations=1;
% Store node locations
node = startPoint:dx:totalLength;



%Random noise bounds
randLow=1;
randHigh=1;


filepath=strcat('stability test v8/rand',num2str(randLow),'/',num2str(numIntervals),'terms');
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


f3prime = @(x) (1/dx - 4/dx^2.*(x) + 3/dx^3.*(x)^2)*dx;
f4prime = @(x) (-2/dx^2.*(x) + 3/dx^3.*(x).^2)*dx;

f3primeprime = @(x) (-4/dx^2 + 6/dx^3.*(x))*dx;
f4primeprime = @(x) (-2/dx^2 + 6/dx^3.*(x))*dx;



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
pGuess = @(x) p1*exp(-gammap1*x)+p2*exp(-gammap2*x)+p3*exp(-gammap3*x);

n1=-1.79577*10^10;
n2=1.29167*10^11;
n3=10^12;
gamman1=100000;
gamman2=10000;
gamman3=1;

%Symbolic form of nGuess
nGuess = @(x) n1*exp(-gamman1*x)+n2*exp(-gamman2*x)+n3*exp(-gamman3*x);

cn1 = Dn - ksi*T*mun/q;
%cn2 = ksi*T*mun/q *diff(pGuess_syms, x)/pGuess_syms;
%cn3 = ksi*T*mun/q*(diff(pGuess_syms, x,2)/pGuess_syms - (diff(pGuess_syms, x)/pGuess_syms)^2) - 1/taun;

cp1 = ksi*T*mup/q - Dp;
%cp2 = -ksi*T*mup/q*diff(nGuess_syms, x)/nGuess_syms;
%cp3 = ksi*T*mup/q*((diff(nGuess_syms, x)/nGuess_syms)^2 - diff(nGuess_syms, x,2)/nGuess_syms) - 1/taup;

nprime = @(x) -cn1.*(n1*gamman1*exp(-gamman1*x) + n2*gamman2*exp(-gamman2*x) + n3*gamman3*exp(-gamman3*x));

Gopn = @(x) cn1.*(n1*(gamman1*gamman1)*exp(-gamman1*x) + n2*(gamman2*gamman2)*exp(-gamman2*x) + n3*(gamman3*gamman3)*exp(-gamman3*x));

fig2 = figure;
plot(node,nGuess(node)), title('nGuess')
saveas(fig2,strcat(filepath,'/nGuess.jpg'))


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
	massmatrix(2*k,2*k) = massmatrix(2*k,2*k) + quad(@(x) f1((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));

	%% phi_right*phi_left
	massmatrix(2*k,2*(k+1)) = massmatrix(2*k,2*(k+1)) + quad(@(x) f1((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));

       
	%% phi_left*phi_right
	massmatrix(2*(k+1),2*(k)) = massmatrix(2*(k+1),2*(k)) + quad(@(x) f2((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));

	%% phi_left*phi_left
	massmatrix(2*(k+1),2*(k+1)) = massmatrix(2*(k+1),2*(k+1)) + quad(@(x) f2((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2)); 

    
    
	%% phi_right*psi_right
	massmatrix(2*k,2*k+1) = massmatrix(2*k,2*k+1) + quad(@(x) f1((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% phi_right*psi_left
	massmatrix(2*k,2*(k+1)+1) = massmatrix(2*k,2*(k+1)+1) + quad(@(x) f1((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% phi_left*psi_right
	massmatrix(2*(k+1),2*k+1) = massmatrix(2*(k+1),2*k+1) + quad(@(x) f2((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% phi_left*psi_left
	massmatrix(2*(k+1),2*(k+1)+1) = massmatrix(2*(k+1),2*(k+1)+1) + quad(@(x) f2((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
    %psi_right*phi_right
    massmatrix(2*k+1,2*k) = massmatrix(2*k+1,2*k) + quad(@(x) f3((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));
    
    %psi_right*phi_left
    massmatrix(2*k+1,2*(k+1))=massmatrix(2*k+1,2*(k+1)) + quad(@(x) f3((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));
    
    
    
    %psi_left*phi_right
    massmatrix(2*(k+1)+1,2*k)=massmatrix(2*(k+1)+1,2*k) + quad(@(x) f4((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));
    
    %psi_left*psi_left
    massmatrix(2*(k+1)+1,2*(k+1))=massmatrix(2*(k+1)+1,2*(k+1)) + quad(@(x) f4((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));  
    
    
    
	%% psi_right*psi_right
	massmatrix(2*k+1,2*k+1) = massmatrix(2*k+1,2*k+1) + quad(@(x) f3((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% psi_right*psi_left
	massmatrix(2*k+1,2*(k+1)+1) = massmatrix(2*k+1,2*(k+1)+1) + quad(@(x) f3((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% psi_left*psi_right
	massmatrix(2*(k+1)+1,2*k+1) = massmatrix(2*(k+1)+1,2*k+1) + quad(@(x) f4((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% psi_left*psi_right
	massmatrix(2*(k+1)+1,2*(k+1)+1) = massmatrix(2*(k+1)+1,2*(k+1)+1) + quad(@(x) f4((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% Phi(k-1) n
	massrhs(2*k) =      massrhs(2*k) +        quad(@(x) f1((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));
	%% Phi(k) n
	massrhs(2*(k+1)) =  massrhs(2*(k+1)) +    quad(@(x) f2((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));
	%% Psi(k-1) n
	massrhs(2*k+1) =    massrhs(2*k+1) +      quad(@(x) f3((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));
	%% Psi(k) n
	massrhs(2*(k+1)+1) = massrhs(2*(k+1)+1) + quad(@(x) f4((x-node(k+1))).*nGuess(x),node(k+1),node(k+2));

end

%Handle 1st interval 'half case'
massmatrix(1,1) = quad(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*(f1((x-node(1))) + alpha*f3((x-node(1)))),node(1),node(1)+dx);
massmatrix(1,2) = quad(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx);
massmatrix(2,1) = quad(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx);
massmatrix(1,3) = quad(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx);
massmatrix(3,1) = quad(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx);

massmatrix(2,2) = massmatrix(2,2) + quad(@(x) f2(x).*f2(x),node(1),node(1)+dx);
massmatrix(2,3) = massmatrix(2,3) + quad(@(x) f2(x).*f4(x),node(1),node(1)+dx);
massmatrix(3,2) = massmatrix(3,2) + quad(@(x) f2(x).*f4(x),node(1),node(1)+dx);
massmatrix(3,3) = massmatrix(3,3) + quad(@(x) f4(x).*f4(x),node(1),node(1)+dx);

massrhs(1) = massrhs(1) + quad(@(x) (f1((x-node(1))) + alpha*f3((x-node(1)))).*nGuess(x),node(1),node(1)+dx);
massrhs(2) = massrhs(2) + quad(@(x) f2((x-node(1))).*nGuess(x),node(1),node(1)+dx);
massrhs(3) = massrhs(3) + quad(@(x) f4((x-node(1))).*nGuess(x),node(1),node(1)+dx);

%Handle last interval 'half case'
last = 2*(numIntervals+1)-2;
massmatrix(last,last) = quad(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*(f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))),node(end)-dx,node(end));
massmatrix(last,last-2) = quad(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-2,last) = quad(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end));

massmatrix(last,last-1) = quad(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-1,last) = quad(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end));

massmatrix(last-2,last-2) = massmatrix(last-2,last-2) + quad(@(x) f1((x-node(end)+dx)).*f1((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-2,last-1) = massmatrix(last-2,last-1) + quad(@(x) f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-1,last-2) = massmatrix(last-1,last-2) + quad(@(x) f1((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end));
massmatrix(last-1,last-1) = massmatrix(last-1,last-1) + quad(@(x) f3((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end));

massrhs(last) = massrhs(last)+ quad(@(x) (f2((x-node(end)+dx)) + beta*f4((x-node(end)+dx))).*nGuess(x),node(end)-dx,node(end));

massrhs(last-2) = massrhs(last-2) + quad(@(x) f1((x-node(end)+dx)).*nGuess(x),node(end)-dx,node(end));
massrhs(last-1) = massrhs(last-1) + quad(@(x) f3((x-node(end)+dx)).*nGuess(x),node(end)-dx,node(end));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform final calcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solve simultaneous equations for basis function coefficients
Coeffsn = massmatrix\massrhs;


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


InitialCnjs = Cnjs;
InitialDnjs = Dnjs;

fig4 = figure;
plot(indexn,Cnjs,node,nGuess(node)),title('nGuess (recomposition)');
saveas(fig4,strcat(filepath,'/nGuess_recomposition.jpg'))


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
	rhs_matrix(2*k) =      rhs_matrix(2*k) +        quad(@(x) f1((x-node(k+1))).*RHS(x),node(k+1),node(k+2));
	%% Phi(k) n
	rhs_matrix(2*(k+1)) =  rhs_matrix(2*(k+1)) +    quad(@(x) f2((x-node(k+1))).*RHS(x),node(k+1),node(k+2));
	%% Psi(k-1) n
	rhs_matrix(2*k+1) =    rhs_matrix(2*k+1) +      quad(@(x) f3((x-node(k+1))).*RHS(x),node(k+1),node(k+2));
	%% Psi(k) n
	rhs_matrix(2*(k+1)+1) = rhs_matrix(2*(k+1)+1) + quad(@(x) f4((x-node(k+1))).*RHS(x),node(k+1),node(k+2));
    
end

rhs_matrix(1) = rhs_matrix(1) + quad(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*RHS(x),node(1),node(1)+dx);
rhs_matrix(2) = rhs_matrix(2) + quad(@(x) f2((x-node(1))).*RHS(x),node(1),node(1)+dx);
rhs_matrix(3) = rhs_matrix(3) + quad(@(x) f4((x-node(1))).*RHS(x),node(1),node(1)+dx);


rhs_matrix(last) = rhs_matrix(last) + quad(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*RHS(x),node(end)-dx,node(end));
rhs_matrix(last-2) = rhs_matrix(last-2) + quad(@(x) f1((x-node(end)+dx)).*RHS(x),node(end)-dx,node(end));
rhs_matrix(last-1) = rhs_matrix(last-1) + quad(@(x) f3((x-node(end)+dx)).*RHS(x),node(end)-dx,node(end));
 
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
	matrix1(2*k,2*k) = matrix1(2*k,2*k) + quad(@(x) f1primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));

	%% phi_right*phi_left
	matrix1(2*k,2*(k+1)) = matrix1(2*k,2*(k+1)) + quad(@(x) f1primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));

       
	%% phi_left*phi_right
	matrix1(2*(k+1),2*(k)) = matrix1(2*(k+1),2*(k)) + quad(@(x) f2primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));

	%% phi_left*phi_left
	matrix1(2*(k+1),2*(k+1)) = matrix1(2*(k+1),2*(k+1)) + quad(@(x) f2primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2)); 

    
    
	%% phi_right*psi_right
	matrix1(2*k,2*k+1) = matrix1(2*k,2*k+1) + quad(@(x) f1primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));%

	%% phi_right*psi_left
	matrix1(2*k,2*(k+1)+1) = matrix1(2*k,2*(k+1)+1) + quad(@(x) f1primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% phi_left*psi_right
	matrix1(2*(k+1),2*k+1) = matrix1(2*(k+1),2*k+1) + quad(@(x) f2primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% phi_left*psi_left
	matrix1(2*(k+1),2*(k+1)+1) = matrix1(2*(k+1),2*(k+1)+1) + quad(@(x) f2primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));%

    
    
    %psi_right*phi_right
    matrix1(2*k+1,2*k) = matrix1(2*k+1,2*k) + quad(@(x)  f3primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));%
    
    %psi_right*phi_left
    matrix1(2*k+1,2*(k+1))=matrix1(2*k+1,2*(k+1)) + quad(@(x) f3primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));
    
    
    
    %psi_left*phi_right
    matrix1(2*(k+1)+1,2*k)=matrix1(2*(k+1)+1,2*k) + quad(@(x) f4primeprime((x-node(k+1))).*f1((x-node(k+1))),node(k+1),node(k+2));
    
    %psi_left*psi_left
    matrix1(2*(k+1)+1,2*(k+1))=matrix1(2*(k+1)+1,2*(k+1)) + quad(@(x) f4primeprime((x-node(k+1))).*f2((x-node(k+1))),node(k+1),node(k+2));  %
    
    
    
	%% psi_right*psi_right
	matrix1(2*k+1,2*k+1) = matrix1(2*k+1,2*k+1) + quad(@(x) f3primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% psi_right*psi_left
	matrix1(2*k+1,2*(k+1)+1) = matrix1(2*k+1,2*(k+1)+1) + quad(@(x) f3primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));

    
    
	%% psi_left*psi_right
	matrix1(2*(k+1)+1,2*k+1) = matrix1(2*(k+1)+1,2*k+1) + quad(@(x) f4primeprime((x-node(k+1))).*f3((x-node(k+1))),node(k+1),node(k+2));

	%% psi_left*psi_right
	matrix1(2*(k+1)+1,2*(k+1)+1) = matrix1(2*(k+1)+1,2*(k+1)+1) + quad(@(x) f4primeprime((x-node(k+1))).*f4((x-node(k+1))),node(k+1),node(k+2));
end

%Handle 1st interval 'half case'
matrix1(1,1) = quad(@(x) (f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*(f1((x-node(1))) + alphan*f3((x-node(1)))),node(1),node(1)+dx);
matrix1(1,2) = quad(@(x) (f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f2((x-node(1))),node(1),node(1)+dx);
matrix1(1,3) = quad(@(x) (f1primeprime((x-node(1))) + alphan*f3primeprime((x-node(1)))).*f4((x-node(1))),node(1),node(1)+dx);
matrix1(2,1) = quad(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*f2primeprime((x-node(1))),node(1),node(1)+dx);
matrix1(3,1) = quad(@(x) (f1((x-node(1))) + alphan*f3((x-node(1)))).*f4primeprime((x-node(1))),node(1),node(1)+dx);

matrix1(2,2) = matrix1(2,2) + quad(@(x) f2primeprime(x).*f2(x),node(1),node(1)+dx);
matrix1(2,3) = matrix1(2,3) + quad(@(x) f2primeprime(x).*f4(x),node(1),node(1)+dx);
matrix1(3,2) = matrix1(3,2) + quad(@(x) f2(x).*f4primeprime(x),node(1),node(1)+dx);
matrix1(3,3) = matrix1(3,3) + quad(@(x) f4primeprime(x).*f4(x),node(1),node(1)+dx);

%Handle last interval 'half case'
matrix1(last,last) =   quad(@(x) (f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*(f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))),node(end)-dx,node(end));
matrix1(last,last-2) = quad(@(x) (f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f1((x-node(end)+dx)),node(end)-dx,node(end));
matrix1(last,last-1) = quad(@(x) (f2primeprime((x-node(end)+dx)) + betan*f4primeprime((x-node(end)+dx))).*f3((x-node(end)+dx)),node(end)-dx,node(end));
matrix1(last-2,last) = quad(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f1primeprime((x-node(end)+dx)),node(end)-dx,node(end));
matrix1(last-1,last) = quad(@(x) (f2((x-node(end)+dx)) + betan*f4((x-node(end)+dx))).*f3primeprime((x-node(end)+dx)),node(end)-dx,node(end));

matrix1(last-2,last-2) = matrix1(last-2,last-2) + quad(@(x) f1((x-node(end)+dx)).*f1primeprime((x-node(end)+dx)),node(end)-dx,node(end));
matrix1(last-2,last-1) = matrix1(last-2,last-1) + quad(@(x) f1primeprime((x-node(end)+dx)).*f3((x-node(end)+dx)),node(end)-dx,node(end));
matrix1(last-1,last-2) = matrix1(last-1,last-2) + quad(@(x) f1((x-node(end)+dx)).*f3primeprime((x-node(end)+dx)),node(end)-dx,node(end));
matrix1(last-1,last-1) = matrix1(last-1,last-1) + quad(@(x) f3((x-node(end)+dx)).*f3primeprime((x-node(end)+dx)),node(end)-dx,node(end));   

    Matrix1Finaln = Cn1*matrix1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %End matrix of n''(x) term
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Vestiges from more expansive versions of code
%Output matrices from these sections should empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iterationCount=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fill matrix associated with n'(x) term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Front end coefficient for this matrix.  DOES NOT INCLUDE CONTRIBUTION FROM
%p(x) DISTRIBUTION.
     Cn2 = ksi*T*mun/q;

     matrix2 = zeros(2*(numIntervals+1)-2,2*(numIntervals+1)-2);

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
plot(indexn,Cnjs,indexn,nGuess(indexn));title('n(x) calculated');
saveas(fig4,strcat(filepath,'/n_calculated_iteration',num2str(iterationCount),'.jpg'))

% fig6 = figure;
% plot(indexn,Dnjs,indexn,subs(diff(nGuess,x,1),x,indexn)),title('nprime calculated');
% saveas(fig6,strcat(filepath,'/nprime_calculated_iteration',num2str(iterationCount),'.jpg'))


 save(strcat(filepath,'/','iteration',num2str(iterationCount),'.mat'))


%Displays runtime of m-file
time = cputime-t

fig10 = figure;
plot_together(InitialCnjs-Cnjs,InitialDnjs-Dnjs,totalLength,numIntervals);
saveas(fig10,strcat(filepath,'/plot_together_iteration',num2str(iterationCount),'.jpg'))
%Store variable values to be used later
save(strcat(filepath,'/',num2str(numIntervals),'terms.mat'))