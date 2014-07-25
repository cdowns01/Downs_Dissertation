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
C=1;

%Set a maximum wavelength to look out to and how frequently to sample.
min1=1100;
min2=1100;
min3=1100;
min4=1100;
max1=4000;
max2=4000;
max3=4000;
max4=4000;
decimation = 50;
%Set up output matrix for ideal efficiencies.  Index of matrix corresponds 
%to a wavelngth (in nm).

%Get photon incidence rate (1/(s*m^2)) for each wavelength
ref3(:,1) = ref(:,1);
for index=1:a-1
    if (ref(index,1)>min1)
    ref3(index,2) = ref(index,column)*(ref(index+1,1)-ref(index,1))/q*ref(index,1)/1240;
                                 %convert to W/m^2            %convert to 1/(s*m^2)                            
    end
end

efficiencies = zeros(max1/decimation-min1/decimation,max2/decimation-min2/decimation,max3/decimation-min3/decimation,max4/decimation-min4/decimation);

%Find efficiency of different band placements
for index=min1/decimation:1:max1/decimation-min1/decimation
    for index2 = min2/decimation:1:max2/decimation-min2/decimation
        for index3 = min3/decimation:1:max3/decimation-min3/decimation
           for index4 = min4/decimation:1:max4/decimation-min4/decimation
                efficiencies(index,index2,index3, index4) = solar_efficiency4(index*decimation+min1,index2*decimation+min2,index3*decimation+min3,index4*decimation+min4,ref3,max_power)/C;
           end
        end
    end
    index
end

%Find maximum efficiency in matrix and extract the indices
[max_efficiency wmax] = max(max(max(max(efficiencies))));
[efficiency zmax] = max(max(max(efficiencies(:,:,:,wmax))));
[efficiency2 ymax] = max(max(efficiencies(:,:,zmax,wmax)));
[efficiency3 xmax] = max(efficiencies(:,ymax,zmax,wmax));

xout = xmax*decimation;  %Most recent value:
x_wave = xout+min1
yout = ymax*decimation;  %Most recent value:
y_wave = yout+min2
zout = zmax*decimation; %Most recent value:
z_wave = zout+min3
wout = wmax*decimation; %Most recent value:
w_wave = wout+min4

max_efficiency  %Most recent value:

%by 5s