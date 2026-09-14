function global_matrix=Assembly_Matrix(global_matrix,local_matrix,LM,e)
% This function is taken from the open-source finite element implementation
% developed by Chennakesava Kadapa [5].
% Original source: https://github.com/chennachaos/ArcLengthMethod

%%% Assemble global_matrix
nen=size(LM,2);
for aa=1:nen
    mm=LM(e,aa);
    if mm~=0
        for bb=1:nen
            nn=LM(e,bb);
            if nn~=0
                global_matrix(mm,nn)=global_matrix(mm,nn)+local_matrix(aa,bb);
            end
        end
    end
end