function H = block_hankel(U, m)
    [M, d] = size(U);

   
    H = zeros(M * m, d - m + 1);  %stack column blocks vertically

    for i = 1:m
        H((i-1)*M + (1:M), :) = U(:, i:d - m + i);
    end
end
