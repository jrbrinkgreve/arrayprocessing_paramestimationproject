function H = generate_H(P)


h = ones(1,P);

for i = 1:P
    if( (i > (P/4)) && (i <= (P/2)) )
        h(i) = -h(i);
    elseif ((i > (3*P/4)) && (i <= (P)) )
        h(i) = -h(i);
    end
end


h = h.';
H = [];


for i = 1:2
    H = blkdiag(H, h);
end


end