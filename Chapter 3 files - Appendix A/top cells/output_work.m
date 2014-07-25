function [work] = output_work(Eg,nph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputs the maximum thermodynamic work a photon can do at a given band
%edge in a solar solar and absorbed photon density(equations from Henry, 1980)
%
%To account for concentration C, change all value of nph to nph*C (or
%include in input argument)
%
%work [eV/photon
%Eg [eV]
%nph [1/(s*m^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e=1.6e-19;
n=3.6; %index of refraction
k=8.617e-5; %[eV/T]
T=300;
hbar = 6.582e-16; %[eVs]
c = 3e8; %[m/s]


A = (e*(n^2+1)*Eg^2*k*T/(4*pi()^2*hbar^3*c^2)); %[A/m^2]

%Calculate Voc
eVoc = Eg-k*T*log(A/(e*nph)); %[eV]

%Initial guess for max power point
eVm_last = eVoc;

%Calculate maximum power point iteratively
for index = 1:10
eVm = eVoc-k*T*log(1+eVm_last/(k*T));
eVm_last=eVm;
end

%Calculate output power in eV
work = Eg - k*T*(log(A/(e*nph)) + log(1+eVm/(k*T)) + 1);