function [region_vert region_list] = anatomy_lookup(MNIvert,region,fname_prefix,opts)
%==========================================================================
% Filename: anatomy_lookup.m (function).
%
% Description:  Function for performing a lookup of which vertices are
%               connected to specific anatomically regions. Can be used to
%               lookup which regions are activated or can be used to
%               build priors based on anatomical structures from an
%               anatomical database
%
% Input:        MNIvert: Vertices coordinates [Nv x 3] in MNI space.
%               region: Cell array with names of regions to be looked up.
%               fname_prefix: Prefix name of files being generate or
%                             loaded - three txt-files and one mat-file:
%                       - [fname_prefix '_NamesOfRegions.txt']
%                       - [fname_prefix '_regions.txt']: Every new line
%                         corresponds to a new region with the first number
%                         specifying the region number of interest.
%                       - [fname_prefix '_vertices.txt']: Every new line
%                         corresponds to a new vertex with the first number
%                         specifying the vertex number of interest.
%                       - [fname_prefix '_regionVSvert.mat']
%               opts: Optional parameters
%                   .database:
%                   .database_path: Path to xjView folder with TDdatabase
%                                   in it.
%                   .load: Load 
%
% Output:       region_vert: Logical array keeping the information of which
%                            vertices being in a given region. Rows:
%                            region number, columns: vertex number.
%               region_list: list with regions specifying which vertices
%                            are included in the given regions
%
% History:
%   - Created:  19/07/2011
%   - Modified: 
%
% Author: Carsten Stahlhut
% Copyright (C) DTU Informatics 2011
%==========================================================================
region_list = [];

%% Check for input variables and optional preferences
if nargin<4, opts.dummy = []; end
if (nargin<3 || isempty(fname_prefix)), fname_prefix = 'List'; end
    
if isfield(opts,'load'), flag_load = opts.load; else flag_load=false; end 

if isfield(opts,'database')
    database = opts.database;
else
    if isfield(opts,'database_path')
        database_path = opts.database_path;
    else
        %     database_path = 'C:\cs\Toolboxes\NeuroTools\SPMadds\xjview';
        database_path = 'C:\Users\Sofie\Dropbox\PhD\Toolboxes\defineRegions';
    end
    database_fname = 'TDdatabase';
    database = load(fullfile(database_path,database_fname));
end

MaskMNIlist = database.wholeMaskMNIAll;
if isempty(region)
%     region = database.DB{6}.anatomy;
    region = fieldnames(MaskMNIlist);
%   mnilist = database.DB{6}.mnilist;
end


%% Now start lookup
Nv = size(MNIvert,1);
Nr = length(region);

if flag_load
    load([fname_prefix '_regionVSvert.mat'])
else
    XYZ = MNIvert;
    region_vert = false(Nr,Nv);
    fid = fopen([fname_prefix '_NamesOfRegions.txt'],'w');
    for ir=1:Nr
        fprintf('Region no. %03.0f out of %.0f is being looked up...',ir,Nr);
        xyz_all = getfield(MaskMNIlist, region{ir});

        for i=1:size(xyz_all,1)
            xyz = xyz_all(i,:);
            [xyz_new,iv,dmin] = spm_XYZreg('NearestXYZ',xyz',XYZ');
            r(i) = iv;
        end
        r = unique(r);
        
%         % Check for outliers - note code not finished for correting for
%         outliers - to be done.
%         Nv_r = length(r);
%         dist = zeros(Nv_r,Nv_r);
%         for i=1:Nv_r
%             dist(:,i) = sum(bsxfun(@minus,XYZ(r(:),:),XYZ(r(i),:)).^2,2);
%         end
%         mean_dist = mean(dist,1);
        
        
        region_vert(ir,r) = true;
        region_list(ir).name = region{ir};
        region_list(ir).ivert = r;

        if ir>1
            dlmwrite([fname_prefix '_regions.txt'],[ir r],'-append','delimiter',',','newline', 'pc');
        else
            dlmwrite([fname_prefix '_regions.txt'],[ir r],'delimiter',',','newline', 'pc');
        end

        fprintf(fid,'%s\n',region{ir});
        fprintf('done.\n');
    end
    fclose(fid);


    for iv=1:Nv
        iv_regions = find(region_vert(:,iv)');
        if ir>1
            dlmwrite([fname_prefix '_vertices.txt'],[iv iv_regions],'-append','delimiter',',','newline', 'pc');
        else
            dlmwrite([fname_prefix '_vertices.txt'],[iv iv_regions],'delimiter',',','newline', 'pc');
        end
    end
    
    save([fname_prefix '_regionVSvert.mat'],'region_vert','region_list')
end




% Nice visualization tool - to be used in the future :o)
Q=load('mri')
Q.D = squeeze(Q.D);
Q.D(:,1:60,:) = [];
p1 = patch(isosurface(Q.D, 5),'FaceColor',[1,.75,.65],...
	'EdgeColor','none');
% p1 = patch(isosurface(D, 5),'FaceColor','red',...
% 	'EdgeColor','none');
p2 = patch(isocaps(Q.D, 5),'FaceColor','interp',...
	'EdgeColor','none');
view(3); axis tight; daspect([1,1,.4])
colormap(gray(100))
camlight left; camlight; lighting gouraud
isonormals(Q.D,p1)





return

%% For testing
clear all, clc, close all
path_headmodel = '/media/FreeAgent GoFlex Drive/cs/Work/Data/ForwardModels/Emotiv';

fname = 'emotiv_SPMleadfields.mat'
database_path = '/media/FreeAgent GoFlex Drive/cs/Toolboxes/NeuroTools/SPMadds/xjview';

D = spm_eeg_load(fullfile(path_headmodel,fname))
val = 1;
MNIvert = D.inv{val}.mesh.tess_mni.vert;
MNIface = D.inv{val}.mesh.tess_mni.face;
fname_prefix = sprintf('List%.0f',size(MNIvert,1))

region = anatomy_brodmann(1:47)
% [region_vert region_list] = anatomy_lookup(MNIvert,region,fname_prefix);
% [region_vert region_list] = anatomy_lookup(MNIvert,[],fname_prefix,struct('load',true));
[region_vert region_list] = anatomy_lookup(MNIvert,region,fname_prefix,struct('database_path',database_path));

[Nr Nv] = size(region_vert);
data = (1:Nr)'*ones(1,Nv);
data = (data.*region_vert)';
for i=1:Nr
    data(:,i) = max(data(:,1:i),[],2);
end
cMAP = hsv(Nr);
cMAP = cMAP(randperm(Nr),:);
handles = plot_3Dbrain(MNIvert,MNIface,data,struct('cMAP',cMAP,'fs',2,'thresh',0,'flag_interp',false));
handles = plot_3Dbrain(MNIvert,MNIface,data,struct('cMAP',cMAP,'fs',2,'thresh',0,'flag_interp',true));
