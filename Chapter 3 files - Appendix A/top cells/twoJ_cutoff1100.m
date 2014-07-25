clear all


ref = xlsread('Am1.5_ASTMG173.xls','SMARTS2');
column =3;

%total_irradiance = sum(ref(:,3))
[a,b] = size(ref);
ref2 = zeros(a,1);
ref3 = zeros(a,2);
q=1.602e-19;

%Find total incident power on the cell
for index=1:a-1
    ref2(index) = ref(index,column)*(ref(index+1,1)-ref(index,1));
                                 %convert to W/m^2                                       
end

max_power = sum(ref2);  %~1000 W/m^2


%Get photon incidence rate (1/(s*m^2)) for each wavelength
ref3(:,1) = ref(:,1);
for index=1:a-1
    ref3(index,2) = ref(index,column)*(ref(index+1,1)-ref(index,1))/q*ref(index,1)/1240;
                                 %convert to W/m^2            %convert to 1/(s*m^2)                            
end

%Set a maximum wavelength to look out to and how frequently to sample.
wavelength_max=1100;
decimation = 1;
%Set up output matrix for ideal efficiencies.  Index of matrix corresponds 
%to a wavelngth (in nm).
efficiencies = zeros(wavelength_max/decimation,wavelength_max/decimation);

%Find efficiency of different band placements
for index=1:1:wavelength_max/decimation
    for index2 = index:1:wavelength_max/decimation
           efficiencies(index,index2) = solar_efficiency2(index*decimation,index2*decimation,ref3,max_power);
    end
    index
end

%Find maximum efficiency in matrix and extract the indices
[max_efficiency ymax] = max(max(efficiencies));
[efficiency xmax] = max(efficiencies(:,ymax));
xout = xmax*decimation %Most recent value:717
yout = ymax*decimation %Most recent value:1099
max_efficiency %Most recent value:42.07%

%by 1ss