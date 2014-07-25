
%Collect Cj's in an array for display purposes
Cnjs = zeros(1,numIntervals+1);
Dnjs = zeros(1,numIntervals+1);
Cpjs = zeros(1,numIntervals+1);
Dpjs = zeros(1,numIntervals+1);

Coeffs=x;
for index=2:length(Cnjs)-1
   Cnjs(1,index) = Coeffs(4*index-5);
   Dnjs(1,index) = Coeffs(4*index-4);
   Cpjs(1,index) = Coeffs(4*index-3);
   Dpjs(1,index) = Coeffs(4*index-2);
end
Cnjs(1) = Coeffs(1);
Dnjs(1) = Coeffs(1)*alphan;
Cpjs(1) = Coeffs(2);
Dpjs(1) = Coeffs(2)*alphap;

Cnjs(numIntervals+1)=Coeffs(length(Coeffs)-1);
Dnjs(numIntervals+1)=Coeffs(length(Coeffs)-1)*betan;
Cpjs(numIntervals+1)=Coeffs(length(Coeffs));
Dpjs(numIntervals+1)=Coeffs(length(Coeffs))*betap;

% Cnjs = Cnjs+delCnjs;
% Cpjs = Cpjs+delCpjs;
% Dnjs = Dnjs+delDnjs;
% Dpjs = Dpjs+delDpjs;
 iterationCount=1;

testCnjs(iterationCount)=norm(Cnjs);
testCpjs(iterationCount)=norm(Cpjs);
testDnjs(iterationCount)=norm(Dnjs);
testDpjs(iterationCount)=norm(Dpjs);

fig4 = figure;
plot(indexn,Cnjs,'x',indexn,nGuess(indexn),indexn,nOrig(indexn));title('n(x) calculated'),legend('n(x) - calculated','n(x) - input','n(x) - original');set(gca,'FontSize',16),xlabel('Depth (in cm)'),ylabel('Carrier concentration (in cm^-3)');
saveas(fig4,strcat(filepath,'/n_calculated_iteration',num2str(iterationCount),'.jpg'))

fig6 = figure;
plot(indexn,Dnjs,'x',indexn,nprime(indexn),indexn,nprimeOrig(indexn)),title('nprime calculated'),legend('Dnjs','nprime - error','nprime - original');
saveas(fig6,strcat(filepath,'/nprime_calculated_iteration',num2str(iterationCount),'.jpg'))

% fig5 = figure;
% plot_together(Cnjs,Dnjs,totalLength,numIntervals);title('plot together');
% saveas(fig5,strcat(filepath,'/plot_together_n',num2str(iterationCount),'.jpg'))

fig4 = figure;
plot(indexn,Cpjs,'x',indexn,pGuess(indexn),indexn,pOrig(indexn));title('p(x) calculated'),legend('p(x) - calculated','p(x) - input','p(x) - orignial');set(gca,'FontSize',16),xlabel('Depth (in cm)'),ylabel('Carrier concentration (in cm^-3)');
saveas(fig4,strcat(filepath,'/p_calculated_iteration',num2str(iterationCount),'.jpg'))

fig6 = figure;
plot(indexn,Dpjs,'x',indexn,pprime(indexn),indexn,pprimeOrig(indexn)),title('pprime calculated'),legend('Dpjs','pGuess - error','pGuess - orignial');
saveas(fig6,strcat(filepath,'/pprime_calculated_iteration',num2str(iterationCount),'.jpg'))
% 
% fig5 = figure;
% plot_together(Cpjs,Dpjs,totalLength,numIntervals);title('plot together');
% saveas(fig5,strcat(filepath,'/plot_together_p',num2str(iterationCount),'.jpg'))