function [hfig cfg] = plot_2Dscalp_useFieldtrip(data,elec,Cnames,opts,cfg)
%==========================================================================
% Filename: plot_2Dscalp.m (function).
% 
% Description:  
%
% Input:        data: Scalp data ( channels x 1 )
%               elec: structure with info of the electrodes. From SPM the
%                     info is obtained by elec = D.sensors('eeg').
%               Cnames: Names of the channels (cell-array)
%               opts:
%                   .hfig: Figure handle if a previous figure should be
%                          used.
%                   .flag_Cnames: Show channels names on plot.
%                   .flag_contour: Show contour of a scalp with a nose.
%
% Output:       hfig: Handles to figure
%
% Example:      hfig = plot_2Dscalp_useFieldtrip(data,elec,Cnames,opts,cfg)
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2011
%==========================================================================

if nargin<4 opts = []; end
if isfield(opts,'hfig'), hfig = opts.hfig; else hfig = figure; end
% if isfield(opts,'flag_Cnames'), flag_Cnames = opts.flag_Cnames; else flag_Cnames = true; end
% if isfield(opts,'flag_contour'), flag_contour = opts.flag_contour; else flag_contour = true; end
if isfield(opts,'time'), time = opts.time; else time = 1:size(data,2); end
if isfield(opts,'colorbar'), cbar = opts.colorbar; else cbar = 'no'; end

% Options of topo plots
if nargin<5
    % cfg.elec = D.sensors('eeg');
    cfg.elec = elec;
    cfg.layout = ft_prepare_layout(cfg);       %Layout of the electrodes
    cfg.xparam = 'time';
    cfg.zparam = 'avg';
    cfg.colorbar = cbar;
    cfg.electrodes = 'labels';
end
% Data info
% data = D(:,:,:);
data_eeg.avg = data;
% data_eeg.time = D.time;
data_eeg.time = time;
data_eeg.dimord = 'chan_time';
% data_eeg.label = D.chanlabels;
data_eeg.label = Cnames;

figure(hfig);
ft_topoplotER(cfg, data_eeg);



