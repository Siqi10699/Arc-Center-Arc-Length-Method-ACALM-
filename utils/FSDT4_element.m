function [Ke,fint,info]=FSDT4_element(Xe,frameE,section,qe,options)
% Analytic residual and consistent tangent for a 4-node FSDT shell, with
% four nodes, five DOFs per node
% The following Kinematics are implemented:
%   * total-Lagrangian degenerated shell
%   * Green-Lagrange strain
%   * finite director rotation through the exponential map
%   * stabilized one-point membrane/bending integration
%   * 1x1 selective reduced integration for transverse shear
%   * energy-consistent full-integration hourglass stabilization
%   * no drilling DOF
if ~isfield(options,'computeTangent')
    options.computeTangent = true;
end
if ~isfield(options,'membraneStabilizationFraction')
    options.membraneStabilizationFraction = 0.08;
end
stabilizationFraction = options.membraneStabilizationFraction;
if ~isscalar(stabilizationFraction) || stabilizationFraction<=0.0 ...
        || stabilizationFraction>1.0
    error(['options.membraneStabilizationFraction must be in (0,1]. ', ...
        'Use 1 for the former full-integration element.']);
end

assert(isequal(size(Xe),[4,3]),'Xe must be 4x3.');
assert(isequal(size(qe),[20,1]),'qe must be 20x1.');
assert(isequal(size(frameE.a1),[4,3]),'frameE.a1 must be 4x3.');
assert(isequal(size(frameE.a2),[4,3]),'frameE.a2 must be 4x3.');
assert(isequal(size(frameE.d0),[4,3]),'frameE.d0 must be 4x3.');

u = zeros(4,3);
d = zeros(4,3);
directorFirst  = zeros(4,3,2);
directorSecond = zeros(4,3,2,2);
for a = 1:4
    ia = 5*(a-1);
    u(a,:) = qe(ia+(1:3)).';
    alpha = qe(ia+4);
    beta  = qe(ia+5);
    [da,d1a,d2a] = director_with_derivatives(alpha,beta,frameE.a1(a,:),frameE.a2(a,:),frameE.d0(a,:));
    d(a,:) = da;
    directorFirst(a,:,:) = d1a;
    directorSecond(a,:,:,:) = d2a;
end
x = Xe+u;

fint = zeros(20,1);
Ke = zeros(20,20);
g = 1.0/sqrt(3.0);
gauss2 = [-g,1.0; g,1.0];

for i = 1:2
    for j = 1:2
        xi = gauss2(i,1);
        eta = gauss2(j,1);
        inPlaneWeight = gauss2(i,2)*gauss2(j,2);
        for p = 1:numel(section.plies)
            ply = section.plies(p);
            [zPoints,zWeights] = thickness_gauss_points(ply.zBottom,ply.zTop);
            for k = 1:2
                [fp,Kp] = integration_point_contribution(Xe,x,frameE,d,directorFirst,directorSecond, ...
                    xi,eta,zPoints(k),[1,2,3],ply.Qm,options.computeTangent);
                weight = stabilizationFraction*inPlaneWeight*zWeights(k);
                fint = fint+weight*fp;
                if options.computeTangent
                    Ke = Ke+weight*Kp;
                end
            end
        end
    end
end

reducedFraction = 1.0-stabilizationFraction;
if reducedFraction>0.0
    for p = 1:numel(section.plies)
        ply = section.plies(p);
        [zPoints,zWeights] = thickness_gauss_points(ply.zBottom,ply.zTop);
        for k = 1:2
            [fp,Kp] = integration_point_contribution(Xe,x,frameE,d,directorFirst,directorSecond, 0.0,0.0,zPoints(k),[1,2,3],ply.Qm,options.computeTangent);
            weight = 4.0*reducedFraction*zWeights(k);
            fint = fint+weight*fp;
            if options.computeTangent
                Ke = Ke+weight*Kp;
            end
        end
    end
end
for p = 1:numel(section.plies)
    ply = section.plies(p);
    [zPoints,zWeights] = thickness_gauss_points(ply.zBottom,ply.zTop);
    for k = 1:2
        [fp,Kp] = integration_point_contribution(Xe,x,frameE,d,directorFirst,directorSecond, 0.0,0.0,zPoints(k),[4,5],ply.Qs,options.computeTangent);
        weight = 4.0*zWeights(k);
        fint = fint+weight*fp;
        if options.computeTangent
            Ke = Ke+weight*Kp;
        end
    end
end

if options.computeTangent
    Ke = 0.5*(Ke+Ke.');
    symmetryError = norm(Ke-Ke.','fro')/max(1.0,norm(Ke,'fro'));
else
    Ke = [];
    symmetryError = NaN;
end

edgeLengths = [norm(Xe(2,:)-Xe(1,:)),norm(Xe(3,:)-Xe(2,:)), ...
               norm(Xe(4,:)-Xe(3,:)),norm(Xe(1,:)-Xe(4,:))];
info.characteristicLength = mean(edgeLengths);
info.tangentSymmetryError = symmetryError;
info.dofOrder = '[ux uy uz alpha beta] per node';
info.residualConvention = 'assemble R = lambda*Fext - Fint';
info.derivativeMethod = 'analytic consistent linearisation';
info.integrationScheme = 'S4R-like stabilized reduced integration';
info.membraneStabilizationFraction = stabilizationFraction;
end

function [fp,Kp] = integration_point_contribution(Xe,x,frameE,d,directorFirst,directorSecond,xi,eta,z,active,D,computeTangent)
[N,dNdxi,dNdeta] = shape_Q4(xi,eta);
referenceLayer = Xe+z*frameE.d0;
currentLayer = x+z*d;
G1 = referenceLayer.'*dNdxi;
G2 = referenceLayer.'*dNdeta;
G3 = frameE.d0.'*N;
g1 = currentLayer.'*dNdxi;
g2 = currentLayer.'*dNdeta;
g3 = d.'*N;

J0 = [G1,G2,G3];
J = [g1,g2,g3];
detJ0 = abs(det(J0));
if detJ0<1.0e-14
    error('Degenerate shell Jacobian. Check connectivity and nodal directors.');
end

e3 = cross(G1,G2);
e3 = e3/norm(e3);
e1 = frameE.a1.'*N;
e1 = e1-e3*(e3.'*e1);
e1 = e1/norm(e1);
e2 = cross(e3,e1);
Q = [e1,e2,e3];

Aref = J0\Q;
H = J*Aref;
Elocal = 0.5*(H.'*H-eye(3));
strain = strain_components(Elocal);
dJ = zeros(3,3,20);
for a = 1:4
    ia = 5*(a-1);
    for component = 1:3
        index = ia+component;
        direction = zeros(3,1);
        direction(component) = 1.0;
        dJ(:,1,index) = dNdxi(a)*direction;
        dJ(:,2,index) = dNdeta(a)*direction;
    end
    for rotation = 1:2
        index = ia+3+rotation;
        direction = reshape(directorFirst(a,:,rotation),3,1);
        dJ(:,1,index) = z*dNdxi(a)*direction;
        dJ(:,2,index) = z*dNdeta(a)*direction;
        dJ(:,3,index) = N(a)*direction;
    end
end

dH = zeros(3,3,20);
Bfull = zeros(5,20);
for index = 1:20
    dH(:,:,index) = dJ(:,:,index)*Aref;
    dE = 0.5*(dH(:,:,index).'*H+H.'*dH(:,:,index));
    Bfull(:,index) = strain_components(dE);
end

B = Bfull(active,:);
stress = D*strain(active);
fp = detJ0*(B.'*stress);

if ~computeTangent
    Kp = [];
    return;
end
Kmaterial = B.'*D*B;

S = stress_tensor(active,stress);
dHvector = reshape(dH,9,20);
Kgeometric = dHvector.'*kron(S,eye(3))*dHvector;

HS = H*S;
for a = 1:4
    ia = 5*(a-1);
    for rotationP = 1:2
        indexP = ia+3+rotationP;
        for rotationQ = 1:rotationP
            indexQ = ia+3+rotationQ;
            direction2 = reshape(directorSecond(a,:,rotationP,rotationQ),3,1);
            d2J = [z*dNdxi(a)*direction2, z*dNdeta(a)*direction2, N(a)*direction2];
            d2H = d2J*Aref;
            value = sum(sum(HS.*d2H));
            Kgeometric(indexP,indexQ) = Kgeometric(indexP,indexQ)+value;
            if indexP~=indexQ
                Kgeometric(indexQ,indexP) = Kgeometric(indexQ,indexP)+value;
            end
        end
    end
end

Kp = detJ0*(Kmaterial+Kgeometric);
end

function value = strain_components(E)
value = [E(1,1);E(2,2);2.0*E(1,2);2.0*E(1,3);2.0*E(2,3)];
end

function S = stress_tensor(active,stress)
S = zeros(3,3);
if isequal(active,[1,2,3])
    S(1,1) = stress(1);
    S(2,2) = stress(2);
    S(1,2) = stress(3);
    S(2,1) = stress(3);
elseif isequal(active,[4,5])
    S(1,3) = stress(1);
    S(3,1) = stress(1);
    S(2,3) = stress(2);
    S(3,2) = stress(2);
else
    error('Unsupported active strain-component set.');
end
end

function [d,dFirst,dSecond] = director_with_derivatives(alpha,beta,a1row,a2row,d0row)
a1 = a1row(:);
a2 = a2row(:);
d0 = d0row(:);
T = [a1,a2];
theta = alpha*a1+beta*a2;
s = theta.'*theta;
[c1,c1s,c1ss,c2,c2s,c2ss] = rodrigues_coefficients(s);
S = skew(theta);
Sp = zeros(3,3,2);
Sp(:,:,1) = skew(a1);
Sp(:,:,2) = skew(a2);
S2 = S*S;
d = d0+c1*S*d0+c2*S2*d0;
dFirst = zeros(3,2);
dSecond = zeros(3,2,2);
sFirst = zeros(2,1);
for p = 1:2
    sFirst(p) = 2.0*(theta.'*T(:,p));
    c1p = c1s*sFirst(p);
    c2p = c2s*sFirst(p);
    dS2p = Sp(:,:,p)*S+S*Sp(:,:,p);
    dFirst(:,p) = c1p*S*d0+c1*Sp(:,:,p)*d0 +c2p*S2*d0+c2*dS2p*d0;
end
for p = 1:2
    for q = 1:2
        sSecond = 2.0*(T(:,p).'*T(:,q));
        c1pq = c1ss*sFirst(p)*sFirst(q)+c1s*sSecond;
        c2pq = c2ss*sFirst(p)*sFirst(q)+c2s*sSecond;
        c1p = c1s*sFirst(p);
        c1q = c1s*sFirst(q);
        c2p = c2s*sFirst(p);
        c2q = c2s*sFirst(q);
        dS2p = Sp(:,:,p)*S+S*Sp(:,:,p);
        dS2q = Sp(:,:,q)*S+S*Sp(:,:,q);
        d2S2 = Sp(:,:,p)*Sp(:,:,q)+Sp(:,:,q)*Sp(:,:,p);
        dSecond(:,p,q) = c1pq*S*d0+c1p*Sp(:,:,q)*d0+c1q*Sp(:,:,p)*d0+c2pq*S2*d0+c2p*dS2q*d0+c2q*dS2p*d0+c2*d2S2*d0;
    end
end
d = d.';
dFirst = permute(dFirst,[3,1,2]);
dSecond = permute(dSecond,[4,1,2,3]);
end

function [c1,c1s,c1ss,c2,c2s,c2ss] = rodrigues_coefficients(s)
if abs(s)<1.0e-4
    c1 = 1-s/6+s^2/120-s^3/5040+s^4/362880-s^5/39916800;
    c1s = -1/6+s/60-s^2/1680+s^3/90720-s^4/7983360;
    c1ss = 1/60-s/840+s^2/30240-s^3/1995840;
    c2 = 1/2-s/24+s^2/720-s^3/40320+s^4/3628800-s^5/479001600;
    c2s = -1/24+s/360-s^2/13440+s^3/907200-s^4/95800320;
    c2ss = 1/360-s/6720+s^2/302400-s^3/23950080;
else
    r = sqrt(s);
    sr = sin(r);
    cr = cos(r);
    c1 = sr/r;
    c1s = (r*cr-sr)/(2*r^3);
    c1ss = (-r^2*sr-3*r*cr+3*sr)/(4*r^5);
    c2 = (1-cr)/r^2;
    c2s = (r*sr-2*(1-cr))/(2*r^4);
    c2ss = (r^2*cr-5*r*sr+8*(1-cr))/(4*r^6);
end
end

function [N,dNdxi,dNdeta] = shape_Q4(xi,eta)
N = 0.25*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
dNdxi = 0.25*[-(1-eta);(1-eta);(1+eta);-(1+eta)];
dNdeta = 0.25*[-(1-xi);-(1+xi);(1+xi);(1-xi)];
end

function [z,w] = thickness_gauss_points(zBottom,zTop)
g = 1.0/sqrt(3.0);
mid = 0.5*(zBottom+zTop);
half = 0.5*(zTop-zBottom);
z = [mid-half*g;mid+half*g];
w = [half;half];
end

function S = skew(v)
S = [0.0,-v(3),v(2);v(3),0.0,-v(1);-v(2),v(1),0.0];
end