function [ifseq,irseq] = insertprimers(seq_i)
    if length(seq_i)>100
        insert_forward = seq_i(1:100);
        ind = length(seq_i)-100;
        insert_reverse = seq_i(ind:end);
        ifseq = insert_forward(1:19);
        irseq = insert_reverse(81:end);
        ifseq = forwarddesign(ifseq,insert_forward);
        irseq = reversedesign(irseq,insert_reverse)
    else
        ifseq = seq_i(1:19);
        ind = length(seq_i)-19;
        irseq = seq_i(ind:end);    
    end
;