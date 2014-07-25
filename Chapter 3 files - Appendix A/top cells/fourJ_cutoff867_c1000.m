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
C=1000;

%Get photon incidence rate (1/(s*m^2)) for each wavelength
ref3(:,1) = ref(:,1);
for index=1:a-1
    ref3(index,2) = C*ref(index,column)*(ref(index+1,1)-ref(index,1))/q*ref(index,1)/1240;
                                 %convert to W/m^2            %convert to 1/(s*m^2)                            
end

%Set a maximum wavelength to look out to and how frequently to sample.
min1=400;
min2=550;
min3=600;
min4=650;
max1=650;
max2=800;
max3=870;
max4=870;
decimation = 5;
%Set up output matrix for ideal efficiencies.  Index of matrix corresponds 
%to a wavelngth (in nm).
efficiencies = zeros(max1/decimation-min1/decimation,max2/decimation-min2/decimation,max3/decimation-min3/decimation,max4/decimation-min4/decimation);

%Find efficiency of different band placements
for index=1:1:max1/decimation-min1/decimation
    for index2 = 1:1:max2/decimation-min2/decimation
        for index3 = 1:1:max3/decimation-min3/decimation
           for index4 = 1:1:max4/decimation-min4/decimation
                efficiencies(index,index2,index3,index4) = solar_efficiency4(index*decimation+min1,index2*decimation+min2,index3*decimation+min3,index4*decimation+min4,ref3,max_power)/C;
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

xout = xmax*decimation;  %Most recent value:525
x_wave = xout+min1
yout = ymax*decimation;  %Most recent value:640
y_wave = yout+min2
zout = zmax*decimation; %Most recent value:750
z_wave = zout+min3
wout = wmax*decimation; %Most recent value:867
w_wave = wout+min4

max_efficiency  %Most recent value:47.26%

%by 5s