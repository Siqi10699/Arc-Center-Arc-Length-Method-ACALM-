function LM=build_FSDT4_LM(elemConn)
%Location matrix for four-node, five-DOF shell elements
nelem = size(elemConn,1);
LM = zeros(nelem,20);
for e = 1:nelem
    for a = 1:4
        node = elemConn(e,a);
        LM(e,5*(a-1)+(1:5)) = 5*(node-1)+(1:5);
    end
end