function zz = reversedesign(reverse_seq, reverse_pool)

if reverse_seq(1) == or('G','C')
    zz = reverse_seq;
else
    next = reverse_pool(length(reverse_pool)-length(reverse_seq)+1);
    while next == 'A'|| next == 'T'
        next = reverse_pool(length(reverse_pool)-length(reverse_seq));
        reverse_seq = strcat(next,reverse_seq);
    end
    zz = reverse_seq;
end

    