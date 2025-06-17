function [angles, freqs] = joint_esprit(data, d, Delta, m)


eps = 1e-8;

%svd to get spatial data in Us and time data in Ut
[Us, ~, ~] = svd(data(:, 1:end-1), 'econ');
[Ut, ~, ~] = svd(transpose(data(1:end-1, :)), 'econ');


%signal subspace up to d
Us = Us(:, 1:d);
Ut = Ut(:, 1:d); 


%smoothing in time means splitting the data in some kind of hankel
%structure:

%we do not use this one atm
Ut_smoothed = block_hankel(Ut, m);
[Ut_smoothed_svU, ~, ~] = svd(Ut_smoothed, 'econ');
Ut_smoothed_svU = Ut_smoothed_svU(:, 1:d);


%{
%smoothing in time means splitting the data in some kind of hankel
%structure:
Ut_smoothed = block_hankel(data.', m);
[Ut_smoothed_svU, ~, ~] = svd(Ut_smoothed, 'econ');
%Ut_smoothed_svU = Ut_smoothed_svU(1:d, :).';
Ut_smoothed_svU = Ut_smoothed_svU(:, 1:d);
%}
%hankel operation not writeable as a matrix multiplication?

%split space data
Usx = Us(1:end-1, :);
Usy = Us(2:end, :);


%split time data
Utx = Ut(1:end-1, :);
Uty = Ut(2:end, :);


Ut_smoothed_svUx = Ut_smoothed_svU(1:end-1, :);
Ut_smoothed_svUy = Ut_smoothed_svU(2:end, :);


%get square matrices to pass to joint_diag
psi_s = pinv(Usx) * Usy;
psi_t = pinv(Utx) * Uty;
psi_smoothed_t = pinv(Ut_smoothed_svUx) * Ut_smoothed_svUy;

%perform joint diagonalisation
[V,D] = joint_diag([psi_s, psi_t], eps);


%retrieve estimates of phases between data points for s and t
angle_phases = diag(D(:,1:d));         %diag(V' * psi_s * V);
freq_phases = diag(D(:,d+1:end));      %diag(V' * psi_smoothed_t * V);


%angle estimation from phase difference in space
angle_rad = asin( angle(angle_phases) / (2 * pi * Delta) );  %invert phi formula
angle_rad = real(angle_rad);                           %take real part
angles = sort(rad2deg(angle_rad));         

%freq estimation from phase difference in time
freqs = sort(mod(angle(freq_phases) / (2*pi), 1));



end

%}



%for some reason the smoothing causes the lines to drift away from another.
%further, the result based on hankel(Ut) vs. hankel(transpose(data)) is not
%the same, which i would have expected as the column spaces are roughly the
%same






%}

%---------------------


