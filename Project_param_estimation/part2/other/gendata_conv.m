function [x, H] = gendata_conv(s, P , N, sigma)

%model the channel
H = generate_H(N, P);

%add noise
noise = sigma * (randn(N*P,1) + 1j*randn(N*P,1)) ./ sqrt(2);


x = H*s + noise;



end
