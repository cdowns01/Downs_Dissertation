function [work1, work2, work3] = output_work3(Eg1, Eg2, Eg3, nph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[work1, work2, work3] = output_work3(Eg1, Eg2, Eg3, nph)
%
%Outputs the maximum thermodynamic work a photon can do at a given band
%edge in a solar solar and absorbed photon density(equations from Henry, 1980)
%Note:nph must be the *lowest* photon density incident on any junction in
%the cell (current matching condition).
%
%To account for concentration C, change all value of nph to nph*C (or
%include in input argument)
%
%work1 [eV/photon absorbed by 1st cell]
%work2 [eV/photon absorbed by 2nd cell]
%Eg [eV]
%nph [1/(s*m^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e=1.6e-19;
n=3.6; %index of refraction
k=8.617e-5; %[eV/T]
T=300;
hbar = 6.582e-16; %[eVs]
c = 3e8; %[m/s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cell 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1 = (e*(n^2+1)*Eg1^2*k*T/(4*pi()^2*hbar^3*c^2)); %[A/m^2]

%Calculate Voc
eVoc1 = Eg1-k*T*log(A1/(e*nph)); %[eV]

%Initial guess for max power point
eVm1_last = eVoc1;

%Calculate maximum power point iteratively
for index = 1:999
eVm1 = eVoc1-k*T*log(1+eVm1_last/(k*T));
eVm1_last=eVm1;
end

%Calcualte output power in eV/photon
work1 = Eg1 - k*T*(log(A1/(e*nph)) + log(1+eVm1/(k*T)) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cell 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A2 = (e*(n^2+1)*Eg2^2*k*T/(4*pi()^2*hbar^3*c^2)); %[A/m^2]

%Calculate Voc
eVoc2 = Eg2-k*T*log(A2/(e*nph)); %[eV]

%Initial guess for max power point
eVm2_last = eVoc2;

%Calculate maximum power point iteratively
for index = 1:999
eVm2 = eVoc2-k*T*log(1+eVm2_last/(k*T));
eVm2_last=eVm2;
end

%Calculate output power in eV/photon
work2 = Eg2 - k*T*(log(A2/(e*nph)) + log(1+eVm2/(k*T)) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cell 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A3 = (e*(n^2+1)*Eg3^2*k*T/(4*pi()^2*hbar^3*c^2)); %[A/m^2]

%Calculate Voc
eVoc3 = Eg3-k*T*log(A3/(e*nph)); %[eV]

%Initial guess for max power point
eVm3_last = eVoc3;

%Calculate maximum power point iteratively
for index = 1:99
eVm3 = eVoc3-k*T*log(1+eVm3_last/(k*T));
eVm3_last=eVm3;
end

%Calculate output power in eV/photon
work3 = Eg3 - k*T*(log(A3/(e*nph)) + log(1+eVm3/(k*T)) + 1);








