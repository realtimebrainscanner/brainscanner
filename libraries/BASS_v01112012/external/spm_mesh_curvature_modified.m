function C = spm_mesh_curvature_modified(M)
% Compute a crude approximation of the curvature of a surface mesh
% FORMAT C = spm_mesh_curvature(M)
% M        - a patch structure
%
% C        - curvature vector
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_curvature.m 3135 2009-05-19 14:49:42Z guillaume $

%CS: problems with structure of M in spm_mesh_adjacency and
%spm_mesh_normals ensures lower and
Mcs = get(M);
Mcs.vertices = Mcs.Vertices;
Mcs.faces = Mcs.Faces;


A = spm_mesh_adjacency(Mcs);
A = sparse(1:size(Mcs.vertices,1),1:size(Mcs.vertices,1),1./sum(A,2)) * A;

C = (A-speye(size(A))) * double(Mcs.vertices);
N = spm_mesh_normals(M);
C = sign(sum(N.*C,2)) .* sqrt(sum(C.*C,2));
