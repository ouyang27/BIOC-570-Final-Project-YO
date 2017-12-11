function seq = loadseq(s)
if length(s) < 10
    accession  = s;
    genbank_dat = getgenbank(accession);
    seq = upper(genbank_dat.Sequence);
else
    seq = upper(s);
end