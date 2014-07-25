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

%Concentration Factor
C=100;

%Get photon incidence rate (1/(s*m^2)) for each wavelength
ref3(:,1) = ref(:,1);
for index=1:a-1
    ref3(index,2) = C*ref(index,column)*(ref(index+1,1)-ref(index,1))/q*ref(index,1)/1240;
                                 %convert to W/m^2            %convert to 1/(s*m^2)                            
end

%Set a maximum wavelength to look out to and how frequently to sample.
min=400;
max1=800;
max2=870;
max3=870;
decimation = 2;
%Set up output matrix for ideal efficiencies.  Index of matrix corresponds 
%to a wavelngth (in nm).
efficiencies = zeros(max1/decimation,max2/decimation,max3/decimation);

%Find efficiency of different band placements
for index=min/decimation:1:max1/decimation
    for index2 = index:1:max2/decimation
        for index3 = index2:1:max3/decimation
                efficiencies(index,index2,index3) = solar_efficiency3(index*decimation,index2*decimation,index3*decimation,ref3,max_power)/C;
        end
    end
    index
end

%Find maximum efficiency in matrix and extract the indices
[max_efficiency zmax] = max(max(max(efficiencies)));
[efficiency ymax] = max(max(efficiencies(:,:,zmax)));
[efficiency2 xmax] = max(efficiencies(:,ymax,zmax));

xout = xmax*decimation  %Most recent value:562
yout = ymax*decimation  %Most recent value:710
zout = zmax*decimation  %Most recent value:870
max_efficiency  %Most recent value:44.64%

%by 2s