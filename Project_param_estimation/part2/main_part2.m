clearvars

nruns = 1000;

wiener_mse_storage = zeros(nruns, 1);
zf_mse_storage = zeros(nruns, 1);


for iteration = 1:nruns

N2 = 500;
P = 4;
sigma = 0.5;


[s_seq, sym] = generate_clean_signal(N2);  % uses N2 for signal length
s_stacked = [s_seq(1:N2-1); s_seq(2:N2)];  % create the stacked data


H = generate_H(P);  %generate the filter matrix
noise = sigma * (randn(2*P, N2-1) + 1j * randn(2*P, N2-1)) / sqrt(2);   %noise matrix so that each row is subject to different noise
received = H*s_stacked + noise; %filter the data and add noise





ZF = pinv(H);     %zero forcing filter with pseudoinverse
WR = (inv((H*H' + (sigma^2)*eye(2*P)))*H)' ;     %wiener filter using the equation in the slides


zf_signal_rec =  ZF(1,:) * received;   %apply the filter using first column
w_signal_rec = WR(1,:) * received;   %apply the filter using first column





%mse calculations

wiener_mse = 0;
zf_mse = 0;
 
for i=1:(N2-1)      %calculate mse of zf and wiener with actual symbols 
    wiener_mse = wiener_mse + abs(w_signal_rec(i)-sym(i))^2;
    zf_mse = zf_mse + abs(zf_signal_rec(i)-sym(i))^2;
end


wiener_mse = wiener_mse/(N2-1);   
zf_mse = zf_mse/(N2-1);   


wiener_mse_storage(iteration) = wiener_mse;
zf_mse_storage(iteration) = zf_mse;


end



P
mean_zf_mse = mean(zf_mse_storage)
mean_wiener_mse = mean(wiener_mse_storage)



%plotting
figure;
hold on;

% Plot ZF recovered symbols in green
scatter(real(zf_signal_rec), imag(zf_signal_rec), 5, 'g', 'filled');
scatter(real(w_signal_rec), imag(w_signal_rec), 5, 'b', 'filled');

% Plot original QPSK symbols in red
scatter(real(sym), imag(sym), 15, 'red', 'filled');

xlabel('In-phase (I)');
ylabel('Quadrature (Q)');
title('QPSK Constellation for P=4: ZF vs Wiener');
legend('ZF', 'Wiener', 'QPSK');
axis equal;
grid on;
hold off;

