function primer = gibsonsingle(vector, insert)
    % get the sequence of the vector
    seq_v = loadseq(vector);

    % get the sequence of the insert gene
    seq_i = loadseq(insert);

    %%

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

    check_seq = seq_v(vsite1-10:vsite1);
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

    [ifseq, irseq] = insertprimers(seq_i);

    %% Ask for extra tag/RE insert as overhang

    extra = input('If any extra insert between the insert gene(s) and the vector, enter 1 for Yes, enter 0 for No: ');
    if extra == 1
        N = input('Please enter the number of extra insert(s): ');
        for ii = 1:N
            fprintf('Extra insert %d  \n',ii);
            seq_tag = input('Please enter the sequence: ');
            seqsite = input('Please enter 0 if inserting on N terminus, 1 if on C terminus: ');
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
                    ir_end = irseq(end-2:end);
                    irseq = irseq(1:end-3);
                    irseq = strcat(irseq,seq_tag,ir_end);
                else
                    irseq = strcat(irseq,seq_tag);
                end
            end
        end
    end

    %% Assembly each primer fragment
    primer.vr.seq = reversecomplm(strcat(vr,ifseq));
    primer.vr.gc = gc(vr);
    primer.vf.seq = strcat(irseq,vf);
    primer.vf.gc = gc(vf);
    primer.if1.seq = strcat(vr,ifseq);
    primer.if1.gc = gc(ifseq);
    primer.ir1.seq = reversecomplm(strcat(irseq,vf));
    primer.ir1.gc = gc(irseq);