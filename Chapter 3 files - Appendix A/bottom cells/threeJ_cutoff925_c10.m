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

%Set a maximum wavelength to look out to and how frequently to sample.
cutoff=925;
min1=926;
min2=1400;
min3=2000;
max1=1300;
max2=1700;
max3=2600;
decimation = 2;

%Concentration Factor
C=10;

for index=1:a-1
    if (ref(index,1)>cutoff)
    ref3(index,2) = C*ref(index,column)*(ref(index+1,1)-ref(index,1))/q*ref(index,1)/1240;
                                 %convert to W/m^2            %convert to 1/(s*m^2)                            
    end                                                      
end

%Set up output matrix for ideal efficiencies.  Index of matrix corresponds 
%to a wavelngth (in nm).
efficiencies = zeros(max1/decimation-min1/decimation,max2/decimation-min2/decimation,max3/decimation-min3/decimation);

%Find efficiency of different band placements
for index=1:1:max1/decimation-min1/decimation
    for index2 = 1:1:max2/decimation-min2/decimation
        for index3 = 1:1:max3/decimation-min3/decimation
                efficiencies(index,index2,index3) = solar_efficiency3(index*decimation+min1,index2*decimation+min2,index3*decimation+min3,ref3,max_power)/C;
        end
    end
    index
end

%Find maximum efficiency in matrix and extract the indices
[max_efficiency zmax] = max(max(max(efficiencies)));
[efficiency ymax] = max(max(efficiencies(:,:,zmax)));
[efficiency2 xmax] = max(efficiencies(:,ymax,zmax));


xout = xmax*decimation;  %Most recent value:1178
x_wave = xout+min1
yout = ymax*decimation;  %Most recent value:1600
y_wave = yout+min2
zout = zmax*decimation; %Most recent value:2436
z_wave = zout+min3

max_efficiency  %Most recent value:13.05%

%By 2s