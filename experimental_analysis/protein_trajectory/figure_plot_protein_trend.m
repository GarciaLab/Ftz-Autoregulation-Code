%% Plot protein trend
% After running part1 of the pipeline, we will use this code to combine all 
% the data together

clear % clear all variables in the workspace
close all % close any figure windows that might be open

addpath(genpath('./lib'))
%% Step 0: Initialization

time = -30:1/6:0;
%% Step 1: Add working path and load the data
% First, we will specify the working path

% Make result path
ResultPath = ['./data/'];
%% 
% Then, read all the data *(not corrected with offset yet)*

load([ResultPath 'all_boundary_traces.mat']);
%% Step 2: Calculate for stripe 4
% use threshold to calculate average concentration

label_stripe = (stripe_num == 4);
protein_input_trace_stripe = protein_input_trace(label_stripe,:);
frame_num_end = length(time);

%%
threshold = 0.5E6;
flag_high = protein_input_trace_stripe(:,frame_num_end)>threshold;
flag_low = protein_input_trace_stripe(:,frame_num_end)<=threshold;
%% Calculate average traces

protein_input_trace_stripe_high = mean(protein_input_trace_stripe(flag_high,:),1);
protein_input_trace_stripe_low = mean(protein_input_trace_stripe(flag_low,:),1);
%% Plot figure

frame_num_start = 121; % -15 min
frame_num_end = size(protein_input_trace_stripe,2);


protein_input_start = protein_input_trace_stripe(:,frame_num_start);
protein_input_end = protein_input_trace_stripe(:,frame_num_end);

protein_input_start(protein_input_start<0) = 0;
protein_input_end(protein_input_end<0) = 0;

protein_input_trace_stripe(protein_input_trace_stripe<0) = 0;


fig = figure;
tiledlayout(1,6)
ax1 = nexttile([1 4]);
hold on
plot(time,protein_input_trace_stripe','Color',[0.75 0.75 0.75])
plot(time,protein_input_trace_stripe(flag_high,:)','Color','#DCECCB')
plot(time,protein_input_trace_stripe(flag_low,:)','Color','#F1D4C9')
plot(time,protein_input_trace_stripe_high','Color','#7AA974','LineWidth',4)
plot(time,protein_input_trace_stripe_low','Color','#D56C55','LineWidth',4)
xlim([-20 0])
ylim([-0.4E6, 2.75E6])


ax2 = nexttile;
histogram(gca,protein_input_start,9,'Normalization','probability') % 9,10 is best?
xlim([-0.4E6, 2.75E6])
set(gca,'view',[90 -90])
set(gcf,'Position',[0 0 600 200])
ylim([0 0.45])


ax3 = nexttile;
histogram(gca,protein_input_end,10,'Normalization','probability') % 10 is best?
xlim([-0.4E6, 2.75E6])
set(gca,'view',[90 -90])
ylim([0 0.5])
set(gcf,'Position',[0 0 600 200])