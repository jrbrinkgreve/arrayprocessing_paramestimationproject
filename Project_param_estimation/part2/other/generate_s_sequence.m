function s = generate_s_sequence(N)


s = zeros(N,1);
for i = 1:N
    s(i) = (sign(rand() - 0.5) +    1j * sign(rand() - 0.5)    ) / sqrt(2);
end



end


