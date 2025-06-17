function [X, A, S] = gendata(M,N,Delta,theta,f,SNR)

%{
Dimensions MxN

M: number of antennas
N: number of samples
Delta: == 1/2

theta (1xd vector): [theta_1 ..... theta_d] directions of sources in deg.
    -90 <= theta_i < 90

f (1xd vector): normalized freqs of each source, 0 <= f_i < 1
SNR (scalar): snr per source, signal power / noise power

s_(i,k) = exp(j 2 pi f_i k), and noise scaled to match snr

-----

%we know x(t) = a*s(t) from the narrowband approximation. we can expand
%this to:


%X = AS + N, with the structure:

% [x1(1)...x1(k)...] = [a1 | a2] * [s1(1).....s1(k)...] + [N....]
% [                ]               [s1(1).....s2(k)...]   [     ]

%}




%run the scripts:
S = generate_S(N, f);
A = generate_A(M, Delta, theta);
noise = generate_N(M, N, S, SNR);
X = A*S + noise;




%-------------------------------------- subfunctions below



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



function noise = generate_N(M, N, S, SNR)


    signal_power = mean(rms(S,2).^2);
    snr_linear = 10^(SNR/10);
    noise_power = signal_power / snr_linear;
    noise = sqrt(noise_power/2) * (randn(M,N) + 1j*randn(M,N));






