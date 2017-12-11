% Primer Design for Gibson Assembly
% Final Project
% Yu Ouyang
%Gibson Assembly for single insert
vector = input('Please enter the accession ID or sequence of your vector: ');
insert = input('Please enter the accession ID(s) or sequence(s) of your insert gene(s): ');
%%
if class(insert) == 'char'
    primer = gibsonsingle(vector, insert);
else 
    primer = gibsonmult(vector, insert);
end
%%
fprintf('Results:  \n');
fprintf('Vector Forward Primer:  \n');
disp(primer.vf)
fprintf('Vector Reverse Primer:  \n');
disp(primer.vr)
if length(insert) > 0 
    fprintf('Insert 1 Forward Primer:  \n');
    disp(primer.if1)
    fprintf('Insert 1 Reverse Primer:  \n');
    disp(primer.ir1)
end
if length(insert) > 1 
    fprintf('Insert 2 Forward Primer:  \n');
    disp(primer.if2)
    fprintf('Insert 2 Reverse Primer:  \n');
    disp(primer.ir2)
end
if length(insert) > 2 
    fprintf('Insert 3 Forward Primer:  \n');
    disp(primer.if3)
    fprintf('Insert 3 Reverse Primer:  \n');
    disp(primer.ir3)
end
if length(insert) > 3 
    fprintf('Insert 4 Forward Primer:  \n');
    disp(primer.if4)
    fprintf('Insert 4 Reverse Primer:  \n');
    disp(primer.ir4)
end
if length(insert) > 4 
    fprintf('Insert 5 Forward Primer:  \n');
    disp(primer.if5)
    fprintf('Insert 5 Reverse Primer:  \n');
    disp(primer.ir5)
end





