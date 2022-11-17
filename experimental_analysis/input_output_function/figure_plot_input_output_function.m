%% Plot input-output function
% After running part1 of the pipeline, we will use this code to combine all 
% the data together

clear % clear all variables in the workspace
close all % close any figure windows that might be open

addpath(genpath('./lib'))
%% Load data

load(['data/all_boundary_traces.mat'],'time','spot_output_trace','protein_input_trace','ap_pos_trace','stripe_num','embryo_num');
%% Plot Input-Output Function
% First, we perform moving average if needed

% specify time for moving average (min)
mean_protein = 1/6;
mean_spot = 1/6;

label_stripe4 = (stripe_num == 4);
protein_input_trace_stripe4 = protein_input_trace(label_stripe4,:);
spot_output_trace_stripe4 = spot_output_trace(label_stripe4,:);

protein_input_trace_stripe4_mean = movmean(protein_input_trace_stripe4,mean_protein*6+1,2,'omitnan');
spot_output_trace_stripe4_mean = movmean(spot_output_trace_stripe4,mean_spot*6+1,2,'omitnan');
%% 
% Then, we need to convert the format

for i = 1:length(time)

    protein_stripe4_temp = protein_input_trace_stripe4_mean(:,i);
    spot_stripe4_temp = spot_output_trace_stripe4_mean(:,i);

    processed_data_stripe4(i).protein = protein_stripe4_temp(~isnan(protein_stripe4_temp));
    processed_data_stripe4(i).spot = spot_stripe4_temp(~isnan(protein_stripe4_temp));
end
%% 
% Then, we can plot the input-output function for each time point

cmap_temp = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560]];
plot_frame = [31, 91, 151];
% 
%% Fit input-output function to Hill Equation
% Binning based on quantiles (stripe 4)

len = 11;

time_plot_start = [-25 -20 -15 -10];
time_plot_end = [-20 -15 -10 -5];

color = ["#D56C55";"#738FC1";"#EAC264";"#7AA974"];

fig = figure;
hold on

for i = 1:4
    time_start = time_plot_start(i);
    time_end = time_plot_end(i);

    io_index = find((time>=time_start) & (time<=time_end));

    x4_temp = [];
    y4_temp = [];
    
    for j = io_index
    
        x4_temp = [x4_temp processed_data_stripe4(j).protein'];
        y4_temp = [y4_temp processed_data_stripe4(j).spot'];
    
    end
    
    x4 = x4_temp(x4_temp>0);
    y4 = y4_temp(x4_temp>0);
    
    x4 = x4';
    y4 = y4';
    
    edges4 = quantile(x4,len-1);
    edges4 = [0 edges4];edges4 = [edges4 max(x4)]; % add zero and max to left and right hand side of the array
    [~,~,loc4]=histcounts(x4,edges4);
    
    yplot4 = accumarray(loc4(loc4>0),y4(loc4>0),[len 1])./accumarray(loc4(loc4>0),1,[len 1]);
    err4 = sqrt(accumarray(loc4(loc4>0),(y4(loc4>0)),[len 1],@(x) mean(x.^2)))./sqrt(accumarray(loc4(loc4>0),1,[len 1]));
    xmid4 = 0.5*(edges4(1:end-1)+edges4(2:end));

   errorbar(xmid4, yplot4, err4, '.-k','MarkerSize',30,'MarkerFaceColor','k','MarkerEdgeColor',color(i),'CapSize',10,'Color',color(i),'LineWidth',1.5);

end

xlim([0 2.25E6])
ylim([0 5.5E5])
xlabel('Protein Concentration (AU)')
ylabel('Transcriptional Output (AU)')
pbaspect([5 9 1])
pbaspect([3 2 1])
%% Generate Fit

%x4(x4<0) = 0;
% generate fit
% specify fit type
fitname = 'a*x^n/(K^n+x^n)';
[io_fit4,gof4] = fit(xmid4(:),yplot4(:),fitname,'StartPoint',[1E6,4E5,3],'Lower',[7E5 3E5 1],'Upper',[2E6,7E5,8]);
[io_fit_full4,gof_full4] = fit(x4(:),y4(:),fitname,'StartPoint',[1E6,4E5,4],'Lower',[7E5 3E5 1],'Upper',[2E6,7E5,8]);

xRange = linspace(0, 2.9E6);


% Plot stripe 4
fig = figure;

hold on

scatter(x4,y4,15,[0.8,0.8,0.8],'filled');
errorbar(xmid4, yplot4, err4, '.k','MarkerSize',20,'MarkerEdgeColor',[0.3,0.3,0.3],'CapSize',15,'LineWidth',1);
plot(xRange,io_fit4(xRange),'LineWidth',2.5);

xlim([0 3E6])
ylim([0 8E5])
xlabel('Ftz protein concentration (au)')
ylabel('transcriptional output (au)')
pbaspect([3 2 1])