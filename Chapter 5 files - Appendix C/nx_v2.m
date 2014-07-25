%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the composition of p(x) for interval preceeding node k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = nx_v2(k,node,Cnjs,Dnjs,f1,f2,f3,f4)

%f1=phi_right
%f2=phi_left
%f3=psi_right
%f4=psi_left
%
out=@(x) Cnjs(k).*f1(x-node(k+1)) + Cnjs(k-1).*f2(x-node(k)) + Dnjs(k).*f3(x-node(k+1)) + Dnjs(k-1).*f4(x-node(k));
% 

% numIntervals = 60;
% totalLength=10^-4;
% startPoint=0;
% intervalLength=totalLength/numIntervals;
% x_axis = linspace((k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,1000);
% 
% figure
% subplot(2,2,1),plot(x_axis,Cnjs(k).*f1(x_axis-node(k+1))),title('phi_right')
% subplot(2,2,2),plot(x_axis,Cnjs(k-1).*f2(x_axis-node(k))),title('phi_left')
% subplot(2,2,3),plot(x_axis,Dnjs(k).*f3(x_axis-node(k+1))),title('psi_right')
% subplot(2,2,4),plot(x_axis,Dnjs(k-1).*f4(x_axis-node(k))),title('psi_left')

end