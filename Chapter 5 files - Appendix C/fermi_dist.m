clear all
close all

E=linspace(-4.5,-3.5,101);
T=300;
Ef = linspace(-4,-3.6,5);
k=8.617*10^-5;

f = zeros(length(E),length(Ef));
for index=1:length(E)
    for findex=1:length(Ef)
f(index,findex) = 1./(1+exp((E(index)-Ef(findex))/(k*T)));
    end
end

plot(E,f),set(gca,'FontSize',20),xlabel('Band Energy (in eV)'),ylabel('Probability state is filled')