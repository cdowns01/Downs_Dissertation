function [efficiency] =  solar_efficiency1(band1, ref3, max_power)
%finds the theoretical efficiency of a solar cell with the input bands
%above
%Input bands are wavelengths (in nm)
%
%
q=1.602e-19;
[a,b]=size(ref3);
sum1 = 0;

%%Find the number of photons incident on the top junction
for index = 1:a-1
   if (ref3(index,1)<=band1)
       sum1 = sum1 + ref3(index,2);
   end
end


band1eV = 1240/band1;
total_work_out = output_work(band1eV,sum1)*q*sum1;
                                        %convert eV to J
efficiency = (total_work_out)/max_power;
