function primer = gibsonmult(vector,insert)
    % get the sequence of the vector
    seq_v = loadseq(vector);

    % get the sequences of insert genes (set up to 5)
    
    seq1 = '';
    seq2 = '';
    seq3 = '';
    seq4 = '';
    seq5 = '';
    if length(insert) > 0
        seq = char(insert(1));
        seq1 = loadseq(seq);
    end
    if length(insert) > 1
        seq = char(insert(2));
        seq2 = loadseq(seq);
    end
    if length(insert) > 2
        seq = char(insert(3));
        seq3 = loadseq(seq);
    end
    if length(insert) > 3
        seq = char(insert(4));
        seq4 = loadseq(seq);
    end
    if length(insert) > 4
        seq = char(insert(5));
        seq5 = loadseq(seq);
    end

    %% Find vector cutting site
        vector_site = input('Please enter the sequence(s) or the position(s) of site(s) where to insert in the vector: ');

    if length(vector_site) == 2 % the insert(s) would replace part of the vector
        if class(vector_site(1)) == 'cell' % input 1st site sequence
            ind = strfind(seq_v,vector_site(1))
            vsite1 = ind+length(cell2mat(vector_site(1)))-1;
            while length(vsite1) > 1 % get rid of non-spocific recognition
                disp('Multiple sequences for position 2 are found in the vector.');
                vector_site(1) = input('please enter new sequence for position 1');
                ind = strfind(seq_v,vector_site(1))
                vsite1 = ind+length(cell2mat(vector_site(1)))-1;
            end
        elseif class(vector_site(1)) == 'double' % input 1st site position
            vsite1 = vector_site(1);
        end
        if class(vector_site(2)) == 'cell' % input 2st site sequence
            vsite2 = strfind(seq_v,vector_site(2));
            while length(vsite2) > 1 
                disp('Multiple sequences for position 2 are found in the vector.');
                vector_site(2) = input('please enter new sequence for position 2');
                vsite2 = strfind(seq_v,vector_site(2));
            end
            vsite2 = vsite2;
        elseif class(vector_site(2)) == 'double' % input 2st site position  
            vsite2 = vector_site(2);
        end
        % need to find which part to be replaced. Mostly, the shorter sequence
        % in the loop would be replaced
        % Default setting: vsite1 for vector reverse primer, vsite2 for vector
        % forward primer
        excision = vsite2 - vsite1;
        para = length(seq_v)./2;
        if excision > para % need to exchange the site if the excision part is outside
            store_vsite1 = vsite1 - length(cell2mat(vector_site(1))) + 1;
            if class(vector_site(2)) == 'char'
                vsite1 = vsite2 + length(cell2mat(vector_site(2))) - 1;
            else
                vsite1 = vsite2;
            end
            vsite2 = store_vsite1;

        elseif class(vector_site(1)) == 'char' 
            vsite1 = vsite1 + length(cell2mat(vector_site(1))) - 1;
            % adjust vsite1 position to include the sequence for position 1 in
            % reverse primer. 
        end

    else % only have one cutting site to insert
        if class(vector_site) == 'char' % input 1st site sequence
            vsite1 = strfind(seq_v,vector_site);
            vsite2 = vsite1;
            vsite1 = vsite1 + length(vector_site) - 1
        else
            vsite1 = vector_site;
            vsite2 = vector_site;
        end
    end


    %% Check if any ATG seq around the vector reverse primer seq, which would affect open reading fragment
    if vsite1>10
        check_seq = seq_v(vsite1-10:vsite1);
    else
        check_seq = strcat(seq_v(end-10:end),seq_v(1:vsite1));
    end
    checkpoint = strfind(check_seq,'ATG');
    if length(checkpoint) > 0
        a = length(check_seq) - checkpoint + 1;
        b = mod(a,3);
        if mod(a,3) == 0 
            vsite1 = vsite1;
        else
            vsite1 = vsite1 + 3 - b
        end
    end


    %% take primer sequences for the vector
    % vector reverse primer vector part, taking 20 bp as default
    L = length(seq_v);
    if vsite1 < 100 % circular structure of vector
        vector_reverse = strcat(seq_v(L-100:end),seq_v(1:vsite1));
    else
        vector_reverse = seq_v(vsite1-100:vsite1);
    end

    vr = vector_reverse(length(vector_reverse)-19:end); %take 20 bp as default
    vr = reversedesign(vr, vector_reverse);

    if vsite2 + 100 > L % circular structure of vector
        vector_forward = strcat(seq_v(vsite2:end),seq_v(1:100));
    else
        vector_forward = seq_v(vsite2:vsite2+100);
    end

    vf = vector_forward(1:19); %take 20 bp as default
    vf = forwarddesign(vf, vector_forward);

    %% get insert primer seq
    [ifseq1, irseq1] = insertprimers(seq1);
    if length(seq2) > 0
        [ifseq2, irseq2] = insertprimers(seq2);
    end
    if length(seq3) > 0
        [ifseq3, irseq3] = insertprimers(seq3);
    end
    if length(seq4) > 0
        [ifseq4, irseq4] = insertprimers(seq4);
    end
    if length(seq5) > 0
        [ifseq5, irseq5] = insertprimers(seq5);
    end
        
    %% Ask for extra tag/RE insert as overhang

    extra = input('If any extra insert between the insert gene(s) and the vector, enter 1 for Yes, enter 0 for No: ');
    if extra == 1
        N = input('Please enter the number of extra insert(s): ');
        for ii = 1:N
            fprintf('Extra insert %d  \n',ii);
            seq_tag = input('Please enter the sequence: ');
            seqnum = input('Please enter which insert gene the extra insert should attach to: ');
            seqsite = input('Please enter 0 if inserting on N terminus, 1 if on C terminus: ');
            if seqnum == 1
                if seqsite == 0
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                        if info == 0
                            vr = strcat(vr,'ATG',seq_tag);
                        else
                            vr = strcat(vr,seq_tag)
                        end
                else
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                    if info == 0
                        ir_end = irseq1(end-2:end);
                        irseq1 = irseq1(1:end-3);
                        irseq1 = strcat(irseq1,seq_tag,ir1_end);
                    else
                        irseq1 = strcat(irseq1,seq_tag);
                    end
                end
            elseif seqnum == 2
                if seqsite == 0
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                        if info == 0
                            irseq1 = strcat(irseq1,'ATG',seq_tag);
                        else
                            irseq1 = strcat(irseq1,seq_tag)
                        end
                else
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                    if info == 0
                        ir2_end = irseq2(end-2:end);
                        irseq2 = irseq2(1:end-3);
                        irseq2 = strcat(irseq2,seq_tag,ir2_end);
                    else
                        irseq2 = strcat(irseq2,seq_tag);
                    end
                end
            elseif seqnum == 3
                if seqsite == 0
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                        if info == 0
                            irseq2 = strcat(irseq2,'ATG',seq_tag);
                        else
                            irseq2 = strcat(irseq2,seq_tag)
                        end
                else
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                    if info == 0
                        ir3_end = irseq3(end-2:end);
                        irseq3 = irseq3(1:end-3);
                        irseq3 = strcat(irseq3,seq_tag,ir3_end);
                    else
                        irseq3 = strcat(irseq3,seq_tag);
                    end
                end
            elseif seqnum == 4
                if seqsite == 0
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                        if info == 0
                            irseq3 = strcat(irseq3,'ATG',seq_tag);
                        else
                            irseq3 = strcat(irseq3,seq_tag)
                        end
                else
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                    if info == 0
                        ir4_end = irseq4(end-2:end);
                        irseq4 = irseq4(1:end-3);
                        irseq4 = strcat(irseq4,seq_tag,ir4_end);
                    else
                        irseq4 = strcat(irseq4,seq_tag);
                    end
                end
            elseif seqnum == 5
                if seqsite == 0
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                        if info == 0
                            irseq4 = strcat(irseq4,'ATG',seq_tag);
                        else
                            irseq4 = strcat(irseq4,seq_tag)
                        end
                else
                    info = input('Please define the extra insert sequence as an extra restriction site (1) or tag sequence (0): '); 
                    if info == 0
                        ir5_end = irseq5(end-2:end);
                        irseq5 = irseq5(1:end-3);
                        irseq5 = strcat(irseq5,seq_tag,ir5_end);
                    else
                        irseq5 = strcat(irseq5,seq_tag);
                    end
                end
            end
        end
    end

    %% Assembly each primer fragment
    primer.vr.seq = reversecomplm(strcat(vr,ifseq1));
    primer.vr.gc = gc(vr);
    primer.if1.seq = strcat(vr,ifseq1);
    primer.if1.gc = gc(ifseq1);
    if length(seq2) > 0
        primer.ir1.seq = reversecomplm(strcat(irseq1,ifseq2));
        primer.ir1.gc = gc(irseq1);
        primer.if2.seq = strcat(irseq1,ifseq2);
        primer.if2.gc = gc(ifseq2);
        if length(seq3) > 0
            primer.ir2.seq = reversecomplm(strcat(irseq2,ifseq3));
            primer.ir2.gc = gc(irseq2);
            primer.if3.seq = strcat(irseq2,ifseq3);
            primer.if3.gc = gc(ifseq3);
            if length(seq4) > 0
                primer.ir3.seq = reversecomplm(strcat(irseq3,ifseq4));
                primer.ir3.gc = gc(irseq3);
                primer.if4.seq = strcat(irseq3,ifseq4);
                primer.if4.gc = gc(ifseq4);
                if length(seq5) > 0
                    primer.ir4.seq = reversecomplm(strcat(irseq4,ifseq5));
                    primer.ir4.gc = gc(irseq4);
                    primer.if5.seq = strcat(irseq4,ifseq5);
                    primer.if5.gc = gc(ifseq5);
                    primer.ir5.seq = reversecomplm(strcat(irseq5,vf));
                    primer.ir5.gc = gc(irseq5);
                    primer.vf.seq = strcat(irseq5,vf);
                    primer.vf.gc = gc(vf);
                else
                    primer.ir4.seq = reversecomplm(strcat(irseq4,vf));
                    primer.ir4.gc = gc(irseq4);
                    primer.vf.seq = strcat(irseq4,vf);
                    primer.vf.gc = gc(vf);
                end
            else
                primer.ir3.seq = reversecomplm(strcat(irseq3,vf));
                primer.ir3.gc = gc(irseq3);
                primer.vf.seq = strcat(irseq3,vf);
                primer.vf.gc = gc(vf);
            end
        else
            primer.ir2.seq = reversecomplm(strcat(irseq2,vf));
            primer.ir2.gc = gc(irseq2);
            primer.vf.seq = strcat(irseq2,vf);
            primer.vf.gc = gc(vf);
        end
    else
        primer.ir1.seq = reversecomplm(strcat(irseq1,vf));
        primer.ir1.gc = gc(irseq1);
        primer.vf.seq = strcat(irseq1,vf);
        primer.vf.gc = gc(vf);
    end
       
    


