function theta = esprit(data, d, Delta)


%svd to get spatial data in Uz
[Uz, ~, ~] = svd(data, 'econ');


%signal subspace up to d
Uz = Uz(:, 1:d);  

%split
Ux = Uz(1:end-1, :);
Uy = Uz(2:end, :);

%pinv and eig
psi = pinv(Ux) * Uy;
lambda = eig(psi);   %get complex exponentials



%angle from psi
angle_rad = asin( angle(lambda) / (2 * pi * Delta) );  %invert phi formula
angle_rad = real(angle_rad);                           %take real part
theta = rad2deg(angle_rad);               



end
