function f = espritfreq(data, d)

%svd column are TIME axis now
data = transpose(data);
[Uz, ~, ~] = svd(data, 'econ');


%signal subspace up to d
Uz = Uz(:, 1:d);  

%split
Ux = Uz(1:end-1, :);
Uy = Uz(2:end, :);

%pinv and eig
psi = pinv(Ux) * Uy;
lambda = eig(psi);   %get complex exponentials


%get freqs from phases        
f = mod(angle(lambda) / (2*pi), 1);  



end

