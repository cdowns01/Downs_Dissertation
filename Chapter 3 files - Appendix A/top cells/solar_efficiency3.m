function [efficiency, ref3, Iout1, Iout2, Iout3] =  solar_efficiency3(band1, band2, band3, ref3, max_power)
%finds the theoretical efficiency of a solar cell with the input bands
%above
%Input bands are wavelengths (in nm)
%
%
q=1.602e-19;

%total_irradiance = sum(ref(:,3))
[a,b] = size(ref3);

sum1 = 0;
for index = 1:a-1
   if (ref3(index,1)<=band1)
       sum1 = sum1 + ref3(index,2);
       ref3(index,2) = 0;
   end
end


sum2 = 0;
for index = 1:a-1
   if (ref3(index,1)<band2 && ref3(index,1)>band1)
       sum2 = sum2 + ref3(index,2);
       ref3(index,2) = 1;
   end
end

sum3 = 0;
for index = 1:a-1
   if (ref3(index,1)<band3 && ref3(index,1)>band2 && ref3(index,1)>band1)
       sum3 = sum3 + ref3(index,2);
       ref3(index,2) = 2;
   end
end

if (sum1<sum2 && sum1<sum3)
   sum_max = sum1;
elseif (sum2<sum3)
   sum_max = sum2;
else
   sum_max = sum3;
end

band1eV = 1240/band1;
band2eV = 1240/band2;
band3eV = 1240/band3;
[cell1_power cell2_power cell3_power] = output_work3(band1eV, band2eV, band3eV, sum_max);

cell1_out = (cell1_power)*q*sum_max;
cell2_out = (cell2_power)*q*sum_max;
cell3_out = (cell3_power)*q*sum_max;

efficiency = (cell1_out+cell2_out+cell3_out)/max_power;
Iout1 = sum1
Iout2 = sum2
Iout3 = sum3
