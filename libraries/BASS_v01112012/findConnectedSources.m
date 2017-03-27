function [Spatch, handlesNew] = findConnectedSources(handles,opts,plotOps)
%==========================================================================
% Filename: findConnectedSources.m (function).
% 
% Description:  
%
% Input:        handles: Object handle
%                   .hfig: Figure handle needs to be specified
%                   .face: Faces in the patch object
%                   .vert: Vertices in the patch object
%               opts:
%                   .vi: Indices on vertices which should be expanded to
%                        patches
%                   .patchSize: Size of desired patch. This number
%                   corresponds to how far from the vi is wanted, i.e.
%                   number of neighbour layers from vi. E.g. patchSize=1
%                   will just give vi, whereas patchSize=2 will give the
%                   layer of the closest surranding neighbors
%                   number of 
%
% Output:       handles: Handles to patch
%
% History:
%   - Created:  28/11/2009
%   - Modified: 23/12/2010
%
% Special remarks: 
%
% Author (C) Carsten Stahlhut, DTU Informatics 2009
%==========================================================================

%
%   handles.hfig: a figure handle needs to be specified
%
% Note this file has been modified 22/12/2010 original file
% findConnectedSources_backup22122010
%
if isfield(opts,'vi')
    vi=opts.vi; 
else
    if isfield(opts,'nSelect'), nSelect = opts.nSelect; else nSelect = 1; end
    vi = SelectVertexPoint(handles.vert,handles.face,nSelect);
end
if isfield(opts,'patchSize'), patchSize = opts.patchSize; else patchSize=3; end

cTess.faces=handles.face;
cTess.vertices=handles.vert;
cTess.vertconn=vertices_connectivity(cTess,0);
% nb=find_vertex_neighbours([],cTess.faces,indx)
% [nb] = find_vertex_neighbours(pnt, tri, indx) %fieldtrip function

Nd = size(handles.vert,1);
Spatch = sparse(Nd,Nd);
for i=1:length(vi)
    Spatch(vi(i),vi(i)) = 1;
    iverts = vi(i);
    for j=2:patchSize
        newverts = expand_patch(iverts, cTess.vertconn);
        iverts = [iverts newverts];
        Spatch(vi(i),newverts) = j;
    end
end

opts_plot.cMAP = lines(patchSize);
opts_plot.thresh = 0;
opts_plot.flag_interp = false;
data = zeros(Nd,1);
for j=patchSize:-1:1
    data(logical(sum(Spatch==j,1)),1) = j;
end
cMAP_brain = [0.5 0.5 0.5; opts_plot.cMAP];
opts_plot.cMAP = cMAP_brain;
opts_plot.cdata = cMAP_brain(data+1,:);

if plotOps == 1
    opts_plot.hfig = handles.hfig;
handlesNew = plot_3Dbrain(handles.vert,handles.face,data,opts_plot);
else 
    handlesNew = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = ConnectedArea(handles,vi,patchSize,cTess)

% if isempty(strfind(pwd,'/'))
%     if ~isempty(strfind(pwd,'P:\'))
%         tempName = D.inv{D.val}.forward.gainmat;
%         [pathstr,name,ext] = fileparts(tempName);
%         pathstr = 'P:\cs\Work\Toolboxes\SPM\Data\multimodal_faces\EEG';
%         D.inv{D.val}.forward.gainmat = [pathstr '\' name ext];
%         
%         tempName = D.inv{D.val}.forward.gainxyz;
%         [pathstr,name,ext] = fileparts(tempName);
%         pathstr = 'P:\cs\Work\Toolboxes\SPM\Data\multimodal_faces\EEG';
%         D.inv{D.val}.forward.gainxyz = [pathstr '\' name ext];        
%     end
% end

% increase patch size:
sc = 1; % swell count init
cp = vi; % current patch init
mag = cell(patchSize); 
mag{sc} = vi;
curve = exp(-((1:patchSize)-1).^2/(2*(patchSize-1)));

ii = patchSize;
while (ii>0)        
    lastpatch{sc}=cp;
    sc=sc+1;
%     ps=patch_swell(cp,cTess.vertconn);  %Appending the next set of
%                                         %adjacent vertices.
    ps = expand_patch(cp, cTess.vertconn);
    cp = [cp(:)' ps];    
    mag{sc} = ps;        
    vc=ones(size(cTess.vertices,1),1)*[0 0 1];
    vc(cp,:)=ones(size(cp,2),1)*[0 1 1];
    set(handles.patch,'FaceVertexCData',vc);
    ii = ii-1;
end
% clear cTess lastpatch ps vc

Nd = size(cTess.vertices,1);
% Nd = size(ATrue,2);
% Nd = 7204;
s = zeros(Nd,1);
for ii=1:length(mag)
    s(mag{ii}) = curve(ii);
end


function ivert_new = expand_patch(iverts, vert_conn)
% PATCH_SWELL: Enlarge a patch by appending the next set of adjacent vertices.
%
% USAGE:  ivert_new = expand_patch(iverts, vert_conn);
%
% INPUT:
%     - iverts : index list of vertex numbers
%     - vert_conn  : sparse matrix of vertex connectivity
% OUTPUT:
%     - ivert_new: row vec list of NEW vertex numbers adjacent the patch



if (size(iverts,1) ~= 1)
  iverts = iverts(:)'; % ensure row vector
end

ivert_new = [];
for i=1:length(iverts)
    ivert_new = [ivert_new vert_conn{iverts(i)}];
end
ivert_new = unique(ivert_new);
% Extract unique set of indices, remove existing vertices
ivert_new = setdiff(ivert_new, iverts);