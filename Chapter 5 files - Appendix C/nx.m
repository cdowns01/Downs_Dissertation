%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the composition of p(x) for interval preceeding node k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = nx(k,Cnjs,Dnjs,phi_right,phi_left,psi_right,psi_left)
% syms x
out= Cnjs(k).*phi_left(k) + Cnjs(k-1).*phi_right(k-1) + Dnjs(k).*psi_left(k) + Dnjs(k-1).*psi_right(k-1);
% 
% phi_left(k)
% numIntervals = 60;
% totalLength=10^-4;
% startPoint=0;
% intervalLength=totalLength/numIntervals;
% x_axis = linspace((k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,1000);
% 
% figure
% subplot(2,2,1),plot(x_axis,subs(Cpjs(k)*phi_left(k),x,x_axis)),title('phi_left')
% subplot(2,2,2),plot(x_axis,subs(Cpjs(k-1)*phi_right(k-1),x,x_axis)),title('phi_right')
% subplot(2,2,3),plot(x_axis,subs(Dpjs(k)*psi_left(k),x,x_axis)),title('psi_left')
% subplot(2,2,4),plot(x_axis,subs(Dpjs(k-1)*psi_right(k-1),x,x_axis)),title('psi_right')
% 
% figure
% plot(x_axis,subs(out,x,x_axis))
end