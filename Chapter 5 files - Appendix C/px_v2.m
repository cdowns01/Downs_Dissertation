%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Returns the composition of p(x) for interval preceeding node k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = px_v2(k,node,Cpjs,Dpjs,f1,f2,f3,f4)
%k = # interval under consideration
%f1=phi_right
%f2=phi_left
%f3=psi_right
%f4=psi_left

out=@(x) Cpjs(k).*f1(x-node(k)) + Cpjs(k+1).*f2(x-node(k)) + Dpjs(k).*f3(x-node(k)) + Dpjs(k+1).*f4(x-node(k));
% 
% Cpjs(k)
% Cpjs(k+1)
% Dpjs(k)
% Dpjs(k+1)
% 
% x_axis = linspace(node(k), node(k+1),100);
% 
% figure
% subplot(2,2,1),plot(x_axis,Cpjs(k).*f1(x_axis-node(k))),title('phi_right')
% subplot(2,2,2),plot(x_axis,Cpjs(k+1).*f2(x_axis-node(k))),title('phi_left')
% subplot(2,2,3),plot(x_axis,Dpjs(k).*f3(x_axis-node(k))),title('psi_right')
% subplot(2,2,4),plot(x_axis,Dpjs(k+1).*f4(x_axis-node(k))),title('psi_left')

end