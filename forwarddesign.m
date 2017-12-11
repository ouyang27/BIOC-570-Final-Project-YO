function zz = forwarddesign(forward_seq,forward_pool)
if forward_seq(length(forward_seq)) == or('G','C')
    zz = forward_seq;
else
    next = forward_pool(length(forward_seq));
    while next == 'A'|| next == 'T'
        next = forward_pool(length(forward_seq)+1);
        forward_seq = strcat(forward_seq,next);
    end
    zz = forward_seq;
end