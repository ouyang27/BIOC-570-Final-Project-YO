function zz = gc(seq)
a = 0;
for ii = 1:length(seq)
    if seq(ii) == 'G' || seq(ii) == 'C'
        a = a+1;
    end
end
zz = a ./ length(seq);