%part 2 of assignment
close all
N = 500;
P = 4;
sigma = 0.5;


%generate a random qpsk sequence
s = generate_s_sequence(N);
[x, H] = gendata_conv(s, P, N, sigma);


X_top = reshape(x(1:end-P), [P, N-1]);
X_bot = reshape(x(P+1:end), [P, N-1]);

X = [X_top; X_bot];


%or can write X = [h 0; 0 h] * hankel(s [2, N-1])




%zero-forcing beamformer: W^H = pinv(H)


%directly on x:

s_hat1 = pinv(H) * x;
rec_error = norm(s_hat1 - s) ./ N


%now putting this in the matrix formulation: (see picture phone 27/05/2025 for derivation)
h = H(1:P, 1).';

W = repmat([h, zeros(1,P)], [N, 1]) ./ (h*h');
s_hat2 = W * X;
rec_error = norm(s_hat2 - s) ./ N

%so i guess this is wrong...


%wiener beamformer: slide 25 lecture spatial processing
W_wiener = inv(H*H' + sigma^2*eye(size(H*H'))) * H;
s_hat_wiener = W_wiener' * x;
rec_error_wiener = norm(s_hat_wiener - s) ./ N





%plotting
figure;
scatter(real(s_hat1)   , imag(s_hat1), '.')
hold on 
scatter(real(s)   , imag(s), '.')
xlim([-1.3 1.3])
ylim([-1.3 1.3])

figure;
scatter(real(s_hat_wiener)   , imag(s_hat_wiener), '.')
hold on 
scatter(real(s)   , imag(s), '.')
xlim([-1.3 1.3])
ylim([-1.3 1.3])






