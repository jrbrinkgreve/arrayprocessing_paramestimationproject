function [spatial_beamformer, temporal_beamformer] =  calc_beamformers(M, Delta, N, angle_estimate, freq_estimate)

%zero-forcing beamformers:

angle_estimate = angle_estimate / 180 * pi; %need to pass angle in rad to function


A_hat = generate_A(M, Delta, angle_estimate);
spatial_beamformer = pinv(A_hat)';



S_hat = generate_S(N, freq_estimate);
temporal_beamformer = pinv(S_hat)';  %need to postmultiply instead of premultiply to get (SS^H) * (SS^H)^-1
                                     %SS^H will be full rank as S is wide
                                     %and full row rank






%-------------- these are just taken from gendata


%generate signal
function S = generate_S(N, f)

S = zeros(length(f), N);

for i = 1:length(f)
    for k = 1:N
        
        S(i,k) = exp(1j*2*pi*f(i)*(k-1));
        
    end
end




function A = generate_A(M, Delta, theta)


%note: a0 = 1, uniform omnidirectional antenna
A = zeros(M, length(theta));
for m = 1:M
    for j = 1:length(theta)

        A(m,j) = exp(1i * (m-1) * 2 * pi * Delta * sin(theta(j)));
%                           ^        "M-1" == m-1 here

    end

end

