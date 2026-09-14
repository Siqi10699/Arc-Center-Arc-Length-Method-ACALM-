function model = hinged_roof_quarter_inputs(materialType,h,nx,nphi)
R = 2540;
L = 254;
theta = 0.1;
nnode = (nx+1)*(nphi+1);
coords = zeros(nnode,3);
nodeId = @(i,j) i*(nphi+1)+j+1;

for i = 0:nx
    x = L*i/nx;
    for j = 0:nphi
        phi = theta*j/nphi;
        node = nodeId(i,j);
        coords(node,:) = [x,R*sin(phi),R*cos(phi)];
    end
end

elemConn = zeros(nx*nphi,4);
e = 0;
for i = 0:nx-1
    for j = 0:nphi-1
        e=e+1;
        elemConn(e,:) = [nodeId(i,j),nodeId(i+1,j),nodeId(i+1,j+1),nodeId(i,j+1)];
    end
end

frames = build_FSDT4_frames(coords,elemConn,[1.0,0.0,0.0]);
switch lower(materialType)
    case {'isotropic','iso'}
        section = FSDT4_section_isotropic(h,3102.75,0.3);
    case {'layup1','0/90/0','layup_0_90_0'}
        section = FSDT4_section_laminate(h,[0;90;0],3300,1100,0.25,660,660,440);
    case {'layup2','90/0/90','layup_90_0_90'}
        section = FSDT4_section_laminate(h,[90;0;90],3300,1100,0.25,660,660,440);
    otherwise
        error('Unknown materialType: %s',materialType);
end
LM = build_FSDT4_LM(elemConn);
neq = 5*nnode;
u = zeros(neq,1);
Fext = zeros(neq,1);
centralNode = nodeId(0,0);
Fext(5*(centralNode-1)+3) = -1000/4;

fixedDofs = [];
for i = 0:nx
    for j = 0:nphi
        node = nodeId(i,j);
        base = 5*(node-1);
        if i==0
            fixedDofs = [fixedDofs,base+1,base+5]; %#ok<AGROW>
        end
        if j==0
            fixedDofs = [fixedDofs,base+2,base+4]; %#ok<AGROW>
        end
        if j==nphi
            fixedDofs = [fixedDofs,base+(1:3)]; %#ok<AGROW>
        end
    end
end
fixedDofs = unique(fixedDofs).';
freeDofs = setdiff((1:neq).',fixedDofs);

model.coords = coords;
model.elemConn = elemConn;
model.frames = frames;
model.sections = {section};
model.LM = LM;
model.u = u;
model.Fext = Fext;
model.fixedDofs = fixedDofs;
model.freeDofs = freeDofs;
model.centralNode = centralNode;
model.dofOrder = '[ux uy uz alpha beta]';
model.fullLoadFromLambda = 'P = 3000*lambda';

model.elemData.frames = frames;
model.elemData.sections = {section};
model.elemData.elementSectionId = ones(size(elemConn,1),1);
model.elemData.elementOptions.computeTangent = true;