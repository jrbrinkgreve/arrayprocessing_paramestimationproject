% Parameters
M = 8; N = 100; d = 2; Delta = 0.5; m = 3;
theta = [30, -10]; % degrees
f = [0.1, 0.25]; % normalized frequency

% Generate synthetic data
angles_rad = deg2rad(theta);
A = exp(1j * 2 * pi * Delta * (0:M-1).' * sin(angles_rad));
S = exp(1j * 2 * pi * f(:) * (0:N-1));
X = A * S + 0.01 * (randn(M,N) + 1j*randn(M,N));

[est_angles, est_freqs] = joint_esprit(X, d, Delta, m);
disp(est_angles)
disp(est_freqs)
