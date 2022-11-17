clc; clear; close all;

fi = ftz_instance();

% load('pooled_embryos_io_all.mat');
% load('pooled_embryos_io_boundary.mat');
% load('pooled_embryos_excluded_io_all.mat');
load('pooled_embryos_excluded_io_boundary.mat');


if ~exist('protein_trace','var')
    if exist('nuclei_fluo_output','var')
         protein_trace = nuclei_fluo_output;
         clear nuclei_fluo_output;
    elseif exist('nuclei_fluo_input','var')
        protein_trace = nuclei_fluo_input;
        clear nuclei_fluo_input;
    end
end

if ~exist('time_point_final','var')
    if exist('time','var')
        time_point_final = time;
        clear time;
    end
end

time_ix = 1:numel(time_point_final);
dRNA_trace = spot_fluo_output;


%% Average traces by binning
% number of nuclei in the dataset
nuclei_num = size(protein_trace,1);
binNum = nuclei_num;
protein_input_avr_bin = protein_trace;

if ~exist('embryoid','var')
    embryoid = ones(size(protein_input_avr_bin,1),1);
end
nembryos = max(embryoid);


%% Smooth the result by spline fitting
protein_input_smooth_bin = nan(binNum,length(time_ix));

for ii = 1:binNum
    cur_val_ix = ~isnan(protein_input_avr_bin(ii,:));
    [curve, ~, ~] = fit(time_point_final(cur_val_ix)', ...
        protein_input_avr_bin(ii,cur_val_ix)', ...
        'smoothingspline','SmoothingParam',0.01);
    protein_input_smooth_bin(ii,cur_val_ix) = curve(time_point_final(cur_val_ix));
end


%% Estimate mRNA dynamics from protein traces
dt = time_point_final(2)-time_point_final(1);

mRNA_est = nan(binNum,length(time_ix));

for ii = 1:binNum
    mRNA_temp = nan(length(time_ix),1);
    protein_temp = protein_input_smooth_bin(ii,:);

    for jj = 1:length(time_ix)-1
        mRNA_temp(jj) = ((protein_temp(jj+1)-protein_temp(jj)) ...
            + fi.gammaP*protein_temp(jj)*dt)/(fi.alpha*dt);
    end
    % Take moving average to smoothen a bit
    % mRNA_mean_temp = movmean(mRNA_temp,13);
    mRNA_mean_temp = mRNA_temp;
    % mRNA value cannot be less than 0 (just in case)
    mRNA_mean_temp(mRNA_mean_temp < 0) = 0;
    mRNA_est(ii,:) = mRNA_mean_temp;
end

mRNA = mRNA_est';
prot = protein_input_smooth_bin';


%% Save
clear fi ii jj cur_val_ix mRNA_est mRNA_temp protein_temp protein_input_smooth_bin

save([batch_name,'_smooth_traces'])