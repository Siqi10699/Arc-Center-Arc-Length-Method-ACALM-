function section = FSDT4_section_laminate(h,anglesDegree,E1,E2,nu12,G12,G13,G23)
% Orthotropic cross-ply/angle-ply shell section
nply = numel(anglesDegree);thicknessFractions = [];shearCorrection = 5/6;
if isempty(thicknessFractions)
    thicknessFractions = ones(nply,1)/nply;
end
if numel(thicknessFractions)~=nply
    error('anglesDegree and thicknessFractions must have equal lengths.');
end
thicknessFractions = thicknessFractions/sum(thicknessFractions);
nu21 = nu12*E2/E1;
den = 1.0-nu12*nu21;
Q11 = E1/den;
Q22 = E2/den;
Q12 = nu12*E2/den;
Q66 = G12;
z = -0.5*h;
plies = repmat(struct('zBottom',0,'zTop',0,'angleDegree',0,'Qm',zeros(3),'Qs',zeros(2)),nply,1);
for p = 1:nply
    angle = anglesDegree(p);
    m = cosd(angle);
    n = sind(angle);
    m2=m*m; n2=n*n; m4=m2*m2; n4=n2*n2;

    Qb11 = Q11*m4+2.0*(Q12+2.0*Q66)*m2*n2+Q22*n4;
    Qb22 = Q11*n4+2.0*(Q12+2.0*Q66)*m2*n2+Q22*m4;
    Qb12 = (Q11+Q22-4.0*Q66)*m2*n2+Q12*(m4+n4);
    Qb16 = (Q11-Q12-2.0*Q66)*m*m2*n-(Q22-Q12-2.0*Q66)*m*n*n2;
    Qb26 = (Q11-Q12-2.0*Q66)*m*n*n2-(Q22-Q12-2.0*Q66)*m*m2*n;
    Qb66 = (Q11+Q22-2.0*Q12-2.0*Q66)*m2*n2+Q66*(m4+n4);
    Qm = [Qb11,Qb12,Qb16; Qb12,Qb22,Qb26; Qb16,Qb26,Qb66];

    shearTransform = [m,n; -n,m];
    Qs = shearCorrection*(shearTransform.'*diag([G13,G23])*shearTransform);
    plyThickness = h*thicknessFractions(p);
    plies(p).zBottom = z;
    plies(p).zTop = z+plyThickness;
    plies(p).angleDegree = angle;
    plies(p).Qm = Qm;
    plies(p).Qs = Qs;
    z = z+plyThickness;
end
section.h = h;
section.plies = plies;
section.type = 'laminate';