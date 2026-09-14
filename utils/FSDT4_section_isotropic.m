function section = FSDT4_section_isotropic(h,E,nu)
% Isotropic plane-stress FSDT shell section
G = E/(2.0*(1.0+nu));shearCorrection = 5/6;
Qm = E/(1.0-nu^2)*[1.0,nu,0.0; nu,1.0,0.0;0.0,0.0,0.5*(1.0-nu)];
Qs = shearCorrection*G*eye(2);
ply.zBottom = -0.5*h;
ply.zTop = 0.5*h;
ply.angleDegree = 0.0;
ply.Qm = Qm;
ply.Qs = Qs;
section.h = h;
section.plies = ply;
section.type = 'isotropic';