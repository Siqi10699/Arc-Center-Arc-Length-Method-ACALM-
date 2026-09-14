function frames=build_FSDT4_frames(coords,elemConn,referenceAxis)
%Build continuous nodal directors and tangent frames
nnode = size(coords,1);
normalSum = zeros(nnode,3);
for e = 1:size(elemConn,1)
    n = elemConn(e,:);
    x = coords(n,:);
    areaVector = cross(x(2,:)-x(1,:),x(3,:)-x(1,:)) + cross(x(3,:)-x(1,:),x(4,:)-x(1,:));
    if norm(areaVector)<1.0e-14
        error('Degenerate element %d.',e);
    end
    for a = 1:4
        normalSum(n(a),:) = normalSum(n(a),:)+areaVector;
    end
end

frames.d0 = zeros(nnode,3);
frames.a1 = zeros(nnode,3);
frames.a2 = zeros(nnode,3);
for node = 1:nnode
    d0 = normalSum(node,:)/norm(normalSum(node,:));
    if size(referenceAxis,1)==1
        trial = referenceAxis;
    else
        trial = referenceAxis(node,:);
    end
    a1 = trial-dot(trial,d0)*d0;
    if norm(a1)<1.0e-10
        candidates = eye(3);
        [~,index] = min(abs(candidates*d0.'));
        trial = candidates(index,:);
        a1 = trial-dot(trial,d0)*d0;
    end
    a1 = a1/norm(a1);
    a2 = cross(d0,a1);
    a2 = a2/norm(a2);
    frames.d0(node,:) = d0;
    frames.a1(node,:) = a1;
    frames.a2(node,:) = a2;
end