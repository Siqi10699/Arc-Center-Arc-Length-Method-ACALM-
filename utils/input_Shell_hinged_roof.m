function [nelem,coords,elemConn,elemData,LM,neq,assy4r,Q]=input_Shell_hinged_roof(materialType,h,nx,nphi)
model = hinged_roof_quarter_inputs(materialType,h,nx,nphi);
nnode = size(model.coords,1);
nelem = size(model.elemConn,1);
coords = model.coords;
elemConn = model.elemConn;
LM = model.LM;
neq = 5*nnode;
assy4r = model.freeDofs;
Q = model.Fext;

elemData.frames = model.frames;
elemData.sections = model.sections;
elemData.elementSectionId = ones(nelem,1);
elemData.elementOptions.computeTangent = true;