%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the composition of p(x) for interval preceeding node k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = px(k,Cpjs,Dpjs,phi_right,phi_left,psi_right,psi_left)
%syms x
out= Cpjs(k).*phi_left(k) + Cpjs(k-1).*phi_right(k-1) + Dpjs(k).*psi_left(k) + Dpjs(k-1).*psi_right(k-1);
% 
% %phi_left(k)
% numIntervals = 100;
% totalLength=10^-4;
% startPoint=0;
% intervalLength=totalLength/numIntervals;
% x_axis = linspace((k-2)*intervalLength+startPoint, (k-1)*intervalLength+startPoint,100);
% 
% % Cpjs(k)
% % phi_left(k)
% % 
% % Cpjs(k-1)
% % phi_right(k-1)
% % 
% % Dpjs(k)
% % psi_left(k)
% % 
% % Dpjs(k-1)
% % psi_right(k-1)
% % 
% % (k-2)*intervalLength+startPoint
% % (k-1)*intervalLength+startPoint
% % 
% figure
% subplot(2,2,1),plot(x_axis,subs(Cpjs(k).*phi_left(k),x,x_axis)),title('phi_left'),xLim([(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint])
% subplot(2,2,2),plot(x_axis,subs(Cpjs(k-1).*phi_right(k-1),x,x_axis)),title('phi_right'),xLim([(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint])
% subplot(2,2,3),plot(x_axis,subs(Dpjs(k).*psi_left(k),x,x_axis)),title('psi_left'),xLim([(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint])
% subplot(2,2,4),plot(x_axis,subs(Dpjs(k-1).*psi_right(k-1),x,x_axis)),title('psi_right'),xLim([(k-2)*intervalLength+startPoint,(k-1)*intervalLength+startPoint])

% figure
% plot(x_axis,subs(out,x,x_axis))
end