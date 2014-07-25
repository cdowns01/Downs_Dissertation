function [interp_data] = plot_together(Cjs,Djs,total_length,numIntervals)

% Set zero array of size to include all nodes and midpoints between
% nodes
 interp_data = zeros(2*numIntervals+1,1);
 
% Copy in nodal data
 interp_data(1:2:end) = Cjs;

% Compute midpoint data
mesh_width = total_length/numIntervals;
 for i=1:numIntervals
 
   phi_left_mid = 1-3/4+ 2/8;
   phi_right_mid = 3/4-2/8;
   psi_left_mid = (1/2- 2/4+1/8)*mesh_width;
   psi_right_mid = -(1/4-1/8)*mesh_width;
   interp_data(2*i) = Cjs(i)*phi_left_mid + Cjs(i+1)*phi_right_mid + Djs(i)*psi_left_mid + Djs(i+1)*psi_right_mid;
 end
% 
plot(0:mesh_width/2:total_length,interp_data);