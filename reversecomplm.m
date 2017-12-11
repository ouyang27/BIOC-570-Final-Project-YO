function zz = reversecomplm(seq)
rcseq = '';
for ii = 1:length(seq)
    if seq(ii) == 'A'
        rc = 'T';
    elseif seq(ii) == 'T'
        rc = 'A';
    elseif seq(ii) == 'G'
        rc = 'C';
    elseif seq(ii) == 'C'
        rc = 'G';
    elseif seq(ii) == 'N'
        rc = 'N';
    else
        rc = seq(ii);
    end
    rcseq = strcat(rc,rcseq);
end
zz = rcseq;
       