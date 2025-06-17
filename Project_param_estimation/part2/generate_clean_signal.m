function [s_seq, sym] = generate_clean_signal(N2)
    s_qpsk = [1/sqrt(2)+1j*1/sqrt(2), -1/sqrt(2)+1j*1/sqrt(2), 1/sqrt(2)-1j*1/sqrt(2), -1/sqrt(2)-1j*1/sqrt(2)];

    s_seq = [];
    sym = [];

    for i = 1:N2
        r = randi([1 4]);
        sym(i) = s_qpsk(r);     % the clean symbol
        s_seq = [s_seq, sym(i)];
    end
end
