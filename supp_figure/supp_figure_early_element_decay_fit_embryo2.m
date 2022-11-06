%% Ftz early element decay dynamics

clear all
close all
clc
%% Parameters

final_frame = 117;
frame_plot = final_frame; % for quality control
fitRange = 57:117; % range for decay fitting
%% Part 1: Read and display data

load(['./data/early_decay_embryo2/early_decay_embryo2_lin.mat']);
load('./data/early_decay_embryo2/CompiledParticles.mat')

%Time = ElapsedTime;
time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
%% Part 2: Quality check on nuclei tracking


% Quality check and save the nuclei that passed the test
% Requires the nuclei to have continuous trace until lateral movement
X_pass = [];
Y_pass = [];
sch_pass = [];

for i = 1:size(schnitzcells,2)
    index = find(schnitzcells(i).frames == final_frame);
    fluo_temp = max(schnitzcells(i).Fluo(1:index,:),[],2);
    result = sum(isnan(fluo_temp));
    % quality check
    if ~isempty(index) && (length(schnitzcells(i).frames)>=final_frame) ...
            && (schnitzcells(i).frames(final_frame) == final_frame) && (result == 0)
        sch_pass = [sch_pass, i];
        X_pass = [X_pass, schnitzcells(i).cenx(index)];
        Y_pass = [Y_pass, schnitzcells(i).ceny(index)];
    end
end

figure(1)
plot(X_pass,Y_pass,'o')
axis equal
xlim([0 512])
ylim([0 512])
title('Nuclei that passed quality check')
%% Part 3: Compile the schnitzcells together

% Initialize storage
processed_data(1).xcoord = [];
processed_data(1).ycoord = [];

processed_data(1).schnitznum = [];
processed_data(1).NuclearFluo = [];
processed_data(1).SpotFluo = [];
processed_data(1).SpotFluoSum = [];
processed_data(1).SpotFluoInt = [];

% Compile all the nuclei in each frame and assign basic info
for i = 1:length(sch_pass)
    for j = 1:final_frame
        sch_now = sch_pass(i);
        frame_now = schnitzcells(sch_now).frames(j);
        x_coord = schnitzcells(sch_now).cenx(j);
        y_coord = schnitzcells(sch_now).ceny(j);
        fluo = max(schnitzcells(sch_now).Fluo(j,:));
        try
            processed_data(frame_now).xcoord = [processed_data(frame_now).xcoord, x_coord];
            processed_data(frame_now).ycoord = [processed_data(frame_now).ycoord, y_coord];
            processed_data(frame_now).schnitznum = [processed_data(frame_now).schnitznum, sch_now];
            processed_data(frame_now).SpotFluo = [processed_data(frame_now).SpotFluo, 0];
            processed_data(frame_now).NuclearFluo = [processed_data(frame_now).NuclearFluo fluo];  
            processed_data(frame_now).SpotFluoSum = [processed_data(frame_now).SpotFluoSum 0];
            processed_data(frame_now).SpotFluoInt = [processed_data(frame_now).SpotFluoInt 0];
        catch
            processed_data(frame_now).xcoord = x_coord;
            processed_data(frame_now).ycoord = y_coord;
            processed_data(frame_now).schnitznum = sch_now;
            processed_data(frame_now).SpotFluo = 0;
            processed_data(frame_now).NuclearFluo = fluo;
            processed_data(frame_now).SpotFluoSum = 0;
            processed_data(frame_now).SpotFluoInt = 0;
        end
    end
end

% degradation constant for mRNA
gamma_r = 1/8; % 1/min

% Assign particle info to all the nuclei
for i = 1:size(CompiledParticles{1,1},2)
    schnitz_num = CompiledParticles{1,1}(i).schnitz;
    for j = 1:size(CompiledParticles{1,1}(i).Frame,2)
        frame = CompiledParticles{1,1}(i).Frame(j);
        if frame<=final_frame
            num = find(processed_data(frame).schnitznum==schnitz_num);
            processed_data(frame).SpotFluo(num) = CompiledParticles{1,1}(i).Fluo(j);
        end
    end
    for j = 1:size(schnitzcells(schnitz_num).frames,1)
        frame = schnitzcells(schnitz_num).frames(j);
        if frame<=final_frame
            num = find(processed_data(frame).schnitznum==schnitz_num);
            % Update integrated spot information
            num_sum = CompiledParticles{1,1}(i).Frame<=frame;
            int_num = sum(CompiledParticles{1,1}(i).Fluo(num_sum));
            processed_data(frame).SpotFluoSum(num) = processed_data(frame).SpotFluoSum(num)+int_num;
            
            % now update the integration value
            frame_par = find(CompiledParticles{1,1}(i).Frame==frame); % find the particle indice at this frame
            if frame>1
                dt = ElapsedTime(frame)-ElapsedTime(frame-1);
                if ~isempty(frame_par)
                    processed_data(frame).SpotFluoInt(num) = processed_data(frame-1).SpotFluoInt(num)*2^(-dt*gamma_r)+...
                        dt*CompiledParticles{1,1}(i).Fluo(frame_par);
                else
                    processed_data(frame).SpotFluoInt(num) = processed_data(frame-1).SpotFluoInt(num)*2^(-dt*gamma_r);
                end
            end
        end
    end
end

xpos = processed_data(frame_plot).xcoord;
ypos = processed_data(frame_plot).ycoord;

% try to plot and test a bit...
pts = [xpos' ypos'];

fig = figure(2);
[v,c] = voronoin(double(pts));

for i = 1:length(c)
    if all(c{i}~=1)
    x = v(c{i},1);
    y = v(c{i},2);
    %a = processed_data(frame_plot).SpotFluo(i);
    a = processed_data(frame_plot).NuclearFluo(i);
    patch(x,y,a);
    colorbar
    %caxis([0 4E6])
    end
end

axis equal
xlim([0 512])
ylim([0 512])

% Convert the format: Store them into individual nuclei traces
nuclei_fluo_traces = zeros(length(sch_pass),final_frame);
sum_spot_fluo_traces = zeros(length(sch_pass),final_frame);
int_spot_fluo_traces = zeros(length(sch_pass),final_frame);
spot_fluo_traces = zeros(length(sch_pass),final_frame);
x_pos = zeros(length(sch_pass),final_frame);
y_pos = zeros(length(sch_pass),final_frame);

for i = 1:length(sch_pass)
    for j = 1:final_frame
        nuclei_fluo_traces(i,j) = processed_data(j).NuclearFluo(i);
        sum_spot_fluo_traces(i,j) = processed_data(j).SpotFluoSum(i);
        int_spot_fluo_traces(i,j) = processed_data(j).SpotFluoInt(i);
        spot_fluo_traces(i,j) = processed_data(j).SpotFluo(i);
        x_pos(i,j) = processed_data(j).xcoord(i);
        y_pos(i,j) = processed_data(j).ycoord(i);
    end
end
%% Part 4: Process traces

% We would like to combine everything together and interpolate values

% Define time point to interpolate
time_point_final = 0:1/6:28;
x = ElapsedTime(1:final_frame);

nuclei_fluo_ins = zeros(length(sch_pass),length(time_point_final));
spot_fluo_ins= zeros(length(sch_pass),length(time_point_final));

nuclei_fluo_output = zeros(length(sch_pass),length(time_point_final));
spot_fluo_input = zeros(length(sch_pass),length(time_point_final));


for i = 1:length(sch_pass)
    
    spot_temp = spot_fluo_traces(i,:);
    int_spot_temp = int_spot_fluo_traces(i,:);
    nuclei_temp = nuclei_fluo_traces(i,:);

    % Interpolate traces
    spot_interp = interp1(x,spot_temp,time_point_final,'pchip');
    int_spot_interp = interp1(x,int_spot_temp,time_point_final,'pchip');
    nuclei_interp = interp1(x,nuclei_temp,time_point_final,'pchip');

    % pmoving average
    spot_interp_avr = movmean(spot_interp,3);
    nuclei_interp_avr = movmean(nuclei_interp,12);

    spot_fluo_ins(i,:) = spot_interp_avr;
    nuclei_fluo_ins(i,:) = nuclei_interp_avr;

    %spot_fluo_input(i,:) = movmean(spot_interp_avr,15);
    % use integrated mRNA level as input
    spot_fluo_input(i,:) = int_spot_interp;
    nuclei_fluo_output(i,:) = movmean(nuclei_interp_avr,15);

    spot_temp = spot_fluo_input(i,:);
    
end
%% Part 5; Estimate decay rate

% Get nuclei for the ones we care about

% Stripe 1
test_frame = 1;
nuclei_stripe = [];

for i = 1:length(sch_pass)
    if (x_pos(i,test_frame)<300) && (x_pos(i,test_frame)>185)
        nuclei_stripe = [nuclei_stripe i];
    end
end

spot_stripe_sum = zeros(final_frame,1);

% Take average of traces
for i = 1:final_frame
    for j = nuclei_stripe
        spot_stripe_sum(i) = spot_stripe_sum(i) + spot_fluo_traces(j,i);
    end
end

%Take moving average
spot_stripe_sum_avr = movmean(spot_stripe_sum,15);

%% Fit this with exponential decay

myFit = fit(time(fitRange)',spot_stripe_sum(fitRange),'exp1');
fig = figure;
hold on

scatter(time(fitRange),spot_stripe_sum(fitRange),30,[0.8,0.8,0.8],'filled');
plot(-25:30,myFit(-25:30),'LineWidth',2)
xlabel('Time (min)')
ylabel('Mean Fluorescence (AU)')
xlim([-25 10])
ylim([0 6E6])
pbaspect([2 3 1])
%%