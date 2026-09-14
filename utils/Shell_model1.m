function [Klocal,Flocal,FintLocal,elementInfo] = Shell_model1(elemData,elemConn,e,coords,u_i)
nodes = shell_element_nodes(elemConn,e);
Xe = coords(nodes,:);
frameE.a1 = elemData.frames.a1(nodes,:);
frameE.a2 = elemData.frames.a2(nodes,:);
frameE.d0 = elemData.frames.d0(nodes,:);
section = shell_element_section(elemData,e);

elementDofs = zeros(20,1);
for a = 1:4
    elementDofs(5*(a-1)+(1:5)) = 5*(nodes(a)-1)+(1:5);
end
qe = u_i(elementDofs);
if isfield(elemData,'elementOptions')
    options = elemData.elementOptions;
else
    options = struct();
end
[Klocal,FintLocal,elementInfo] = FSDT4_element(Xe,frameE,section,qe,options);

Flocal = -FintLocal;
end

function nodes = shell_element_nodes(elemConn,e)
if size(elemConn,2)==4
    nodes = elemConn(e,1:4);
elseif size(elemConn,2)>=6
    nodes = elemConn(e,end-3:end);
else
    error('Shell elemConn must contain four node labels per element.');
end
nodes = double(nodes);
end

function section = shell_element_section(elemData,e)
sections = elemData.sections;
if isfield(elemData,'elementSectionId')
    sectionId = elemData.elementSectionId(e);
else
    sectionId = 1;
end
if iscell(sections)
    section = sections{sectionId};
elseif numel(sections)==1
    section = sections;
else
    section = sections(sectionId);
end
end