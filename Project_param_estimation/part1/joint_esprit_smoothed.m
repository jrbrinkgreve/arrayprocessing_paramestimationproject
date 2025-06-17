function [angles, freqs] = joint_esprit_smoothed(data, d, Delta, m)

eps = 1e-8;

[M, N] = size(data);

hankelized_data = block_hankel(data, m);


[U_temp, ~ ,~] = svd(hankelized_data, 'econ');
U = U_temp(:,1:d); %d signals truncate


%time // angle estim
% J_xphi: selects the first (m-1) blocks of M rows
%size: [M*(m-1), M*m]
%this selecs the top m-1 and bottom m-1 U-s
Jx_phi = kron([eye(m-1), zeros(m-1,1)], eye(M));
Jy_phi = kron([zeros(m-1,1), eye(m-1)], eye(M));

%space
% size: [(M-1)*m, M*m]
%this selects the top and bottom part of each block stacked U, first and last M-1
%rows.
Jx_theta = kron(eye(m), [eye(M-1), zeros(M-1,1)]);
Jy_theta = kron(eye(m), [zeros(M-1,1), eye(M-1)]);


%time pinvs
U_x_phi = Jx_phi * U;
U_y_phi = Jy_phi * U;


%space pinvs
U_x_theta = Jx_theta * U;
U_y_theta = Jy_theta * U;

%similarity transform stuff
psi_t = pinv(U_x_phi) * U_y_phi;
psi_s = pinv(U_x_theta) * U_y_theta;


%perform joint diagonalisation
[V,D] = joint_diag([psi_s, psi_t], eps);


%retrieve estimates of phases between data points for s and t
angle_phases = diag(D(:,1:d));         %diag(V' * psi_s * V);
freq_phases = diag(D(:,d+1:2*d));      %diag(V' * psi_smoothed_t * V);


%angle estimation from phase difference in space
angle_rad = asin( angle(angle_phases) / (2 * pi * Delta) );  %invert phi formula
angle_rad = real(angle_rad);                           %take real part
angles = sort(rad2deg(angle_rad));         

%freq estimation from phase difference in time
freqs = sort(mod(angle(freq_phases) / (2*pi), 1));



end





