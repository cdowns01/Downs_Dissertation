clear all


ref = xlsread('Am1.5_ASTMG173.xls','SMARTS2');
column =3;

%total_irradiance = sum(ref(:,3))
[a,b] = size(ref);
ref2 = zeros(a,1);
ref3 = zeros(a,2);
q=1.602e-19;

for index=1:a-1
    ref2(index) = ref(index,column)*(ref(index+1,1)-ref(index,1));
                                 %convert to W/m^2                                       
end

max_power = sum(ref2);  %~1000 W/m^2

%Concentration factor
C=1000;

ref3(:,1) = ref(:,1);
%Get photon incidence rate (1/(s*m^2)) for each wavelength
for index=1:a-1
    ref3(index,2) = C*ref(index,column)*(ref(index+1,1)-ref(index,1))/q*ref(index,1)/1240;
                                 %convert to W/m^2            %convert to 1/(s*m^2)                            
end

wavelength_max=867;
decimation = 1;

efficiencies = zeros(1,1500);

for index=1:1:wavelength_max/decimation
   efficiencies(index) = solar_efficiency1(index*decimation,ref3,max_power)/C;
   index
end

[max_efficiency xmax] = max(efficiencies);
xout = xmax*decimation %Most recent value:867
max_efficiency %Most recent value:36.46%