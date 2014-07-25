function [efficiency] =  solar_efficiency2(band1, band2, ref3, max_power)
%finds the theoretical efficiency of a solar cell with the input bands
%above
%Input bands are wavelengths (in nm)
%
%
q=1.602e-19;
[a,b] = size(ref3);


sum1 = 0;
for index = 1:a-1
   if (ref3(index,1)<=band1)
       sum1 = sum1 + ref3(index,2);
   end
end


sum2 = 0;
for k = 1:a-1
   if (ref3(k,1)<band2 && ref3(k,1)>band1)
       sum2 = sum2 + ref3(k,2);
   end
end

if (sum1<sum2)
   sum_max = sum1;
else
   sum_max = sum2;
end

band1eV = 1240/band1;
band2eV = 1240/band2;
[cell1_power cell2_power] = output_work2(band1eV, band2eV, sum_max);

cell1_out = (cell1_power)*q*sum_max;
cell2_out = (cell2_power)*q*sum_max;

efficiency = (cell1_out+cell2_out)/max_power;
