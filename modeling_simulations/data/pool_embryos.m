clc; clear; close all;

folder = 'embryos_for_input_output_function';
batch_name = 'pooled_embryos_io_all';
files = dir(fullfile('.',filesep,folder,filesep,'2*all*.mat'));

% folder = 'embryos_for_input_output_function';
% batch_name = 'pooled_embryos_io_boundary';
% files = dir(fullfile('.',filesep,folder,filesep,'2*boundary*.mat'));

% folder = 'embryos_excluded_from_input_output_function';
% batch_name = 'pooled_embryos_excluded_io_all';
% files = dir(fullfile('.',filesep,folder,filesep,'2*all*.mat'));

% folder = 'embryos_excluded_from_input_output_function';
% batch_name = 'pooled_embryos_excluded_io_boundary';
% files = dir(fullfile('.',filesep,folder,filesep,'2*boundary*.mat'));

files = {files.name};

nembryos = numel(files);

for ii = 1:nembryos
    cur(ii) = load(files{ii});
end

%%
[max_time_pts,max_ix] = max(arrayfun(@(x) length(x.time_point_final),cur));
nnuclei = sum(arrayfun(@(x) numel(x.xcoord),cur));

nuclei_fluo_input = nan(nnuclei,max_time_pts);
spot_fluo_output = nan(nnuclei,max_time_pts);
time_point_final = cur(max_ix).time_point_final;
xcoord = nan(nnuclei,1);
ycoord = nan(nnuclei,1);
embryoid = nan(nnuclei,1);

start_ix_nuc = 1;
for ii = 1:nembryos
    cur_nnuclei = numel(cur(ii).xcoord);
    cur_ixs = start_ix_nuc + (0:cur_nnuclei-1);
    embryoid(cur_ixs) = ii;
    
    % assume all traces end at time 0
    curntpts = numel(cur(ii).time_point_final);
    cur_tixs = max_time_pts + ((-curntpts+1):0);
    
    nuclei_fluo_input(cur_ixs,cur_tixs) = cur(ii).nuclei_fluo_input;
    spot_fluo_output(cur_ixs,cur_tixs) = cur(ii).spot_fluo_output;
    xcoord(cur_ixs) = cur(ii).xcoord;
    ycoord(cur_ixs) = cur(ii).ycoord;
    
    start_ix_nuc = start_ix_nuc + cur_nnuclei;
end

%%
clear cur cur_ixs cur_nnuclei cur_tixs curntpts ii max_ix ...
    max_time_pts files start_ix_nuc

save(batch_name)


% %%
% figure
% for cur_embryoid = 1:nembryos
%     cur_emb_ix = (embryoid == cur_embryoid);
%     x_midpoint = mean([min(xcoord(cur_emb_ix)),max(xcoord(cur_emb_ix))]);
%     stripe_4_ix = cur_emb_ix & (xcoord < x_midpoint);
%     stripe_5_ix = cur_emb_ix & (xcoord > x_midpoint);
%     
%     subplot(1,nembryos,cur_embryoid)
%     scatter(xcoord(stripe_4_ix),ycoord(stripe_4_ix),150,'filled'); hold on;
%     scatter(xcoord(stripe_5_ix),ycoord(stripe_5_ix),150,'filled'); hold off;
%     axis off
% end
% set(gcf,'Position',[60 520 300*nembryos,225]);