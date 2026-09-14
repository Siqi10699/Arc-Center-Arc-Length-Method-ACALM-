function global_vector=Assembly_Vector(global_vector,local_vector,LM,e)
% This function is taken from the open-source finite element implementation
% developed by Chennakesava Kadapa [5].
% Original source: https://github.com/chennachaos/ArcLengthMethod

%%% Assemble F_global
nen=size(LM,2);
for aa=1:nen
    mm=LM(e,aa);
    if mm~=0
        global_vector(mm,1)=global_vector(mm,1)+local_vector(aa,1);
    end
end