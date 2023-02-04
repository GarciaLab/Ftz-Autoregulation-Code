%% Ftz early element expression dynamics

clear all
close all
clc

spot_fluo_full = [];
%% loop through three embryos

for m = 1:7
%% Part 1: Read and display data

switch m
    case 1

        final_frame = 95;
        frame_plot = final_frame; % for quality control
        
        load('./data/late_dynamics_yellow_embryo1/late_dynamics_yellow_embryo1_lin.mat');
        load('./data/late_dynamics_yellow_embryo1/CompiledParticles.mat')
    
        time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
    
        x_pos_max = 335;
        x_pos_min = 200;

    case 2

        final_frame = 90;
        frame_plot = final_frame; % for quality control
        
        load('./data/late_dynamics_yellow_embryo2/late_dynamics_yellow_embryo2_lin.mat');
        load('./data/late_dynamics_yellow_embryo2/CompiledParticles.mat')
    
        time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
    
        x_pos_max = 340;
        x_pos_min = 200;

    case 3

        final_frame = 56;
        frame_plot = final_frame; % for quality control
        
        load('./data/late_dynamics_yellow_embryo3/late_dynamics_yellow_embryo3_lin.mat');
        load('./data/late_dynamics_yellow_embryo3/CompiledParticles.mat')
    
        time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
    
        x_pos_max = 350;
        x_pos_min = 210;

    case 4

        final_frame = 61;
        frame_plot = final_frame; % for quality control
        
        load('./data/late_dynamics_yellow_embryo4/late_dynamics_yellow_embryo4_lin.mat');
        load('./data/late_dynamics_yellow_embryo4/CompiledParticles.mat')
    
        time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
    
        x_pos_max = 310;
        x_pos_min = 180;

    case 5

        final_frame = 80;
        frame_plot = final_frame; % for quality control
        
        load('./data/late_dynamics_yellow_embryo5/late_dynamics_yellow_embryo5_lin.mat');
        load('./data/late_dynamics_yellow_embryo5/CompiledParticles.mat')
    
        time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
    
        x_pos_max = 325;
        x_pos_min = 200;

    case 6

        final_frame = 64;
        frame_plot = final_frame; % for quality control
        
        load('./data/late_dynamics_yellow_embryo6/late_dynamics_yellow_embryo6_lin.mat');
        load('./data/late_dynamics_yellow_embryo6/CompiledParticles.mat')
    
        time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
    
        x_pos_max = 340;
        x_pos_min = 200;

    case 7

        final_frame = 70;
        frame_plot = final_frame; % for quality control
        
        load('./data/late_dynamics_yellow_embryo7/late_dynamics_yellow_embryo7_lin.mat');
        load('./data/late_dynamics_yellow_embryo7/CompiledParticles.mat')
    
        time = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);
    
        x_pos_max = 320;
        x_pos_min = 170;

end


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
            && (schnitzcells(i).frames(final_frame) == final_frame) && (result == 0) && any(fluo_temp)
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

processed_data = [];

% Initialize storage
processed_data(1).xcoord = [];
processed_data(1).ycoord = [];

processed_data(1).schnitznum = [];
processed_data(1).SpotFluo = [];

sch_final_num = 0;

% Compile all the nuclei in each frame and assign basic info
for i = 1:length(sch_pass)

    sch_now = sch_pass(i);
    x_avr = mean(schnitzcells(sch_now).cenx);

    if (x_avr<=x_pos_max) & (x_avr>=x_pos_min)
        sch_final_num = sch_final_num + 1;

        for j = 1:final_frame
            frame_now = schnitzcells(sch_now).frames(j);
            x_coord = schnitzcells(sch_now).cenx(j);
            y_coord = schnitzcells(sch_now).ceny(j);
            fluo = max(schnitzcells(sch_now).Fluo(j,:));
            try
                processed_data(frame_now).xcoord = [processed_data(frame_now).xcoord, x_coord];
                processed_data(frame_now).ycoord = [processed_data(frame_now).ycoord, y_coord];
                processed_data(frame_now).schnitznum = [processed_data(frame_now).schnitznum, sch_now];
                processed_data(frame_now).numSpot = [processed_data(frame_now).numSpot 0];
                processed_data(frame_now).SpotFluoSum =[processed_data(frame_now).SpotFluoSum 0];
            catch
                processed_data(frame_now).xcoord = x_coord;
                processed_data(frame_now).ycoord = y_coord;
                processed_data(frame_now).schnitznum = sch_now;
                processed_data(frame_now).numSpot = 0;
                processed_data(frame_now).SpotFluoSum = 0;
            end
        end
    end
end


% Assign particle info to all the nuclei
for i = 1:size(CompiledParticles{1,1},2)
    schnitz_num = CompiledParticles{1,1}(i).schnitz;
    for j = 1:size(CompiledParticles{1,1}(i).Frame,2)
        frame = CompiledParticles{1,1}(i).Frame(j);
        if frame<=final_frame
            num = find(processed_data(frame).schnitznum==schnitz_num);
            processed_data(frame).numSpot(num) = processed_data(frame).numSpot(num)+1;
            processed_data(frame).SpotFluoSum(num) = processed_data(frame).SpotFluoSum(num)+CompiledParticles{1,1}(i).Fluo(j);
        end
    end
end

xpos = processed_data(frame_plot).xcoord;
ypos = processed_data(frame_plot).ycoord;

axis equal
xlim([0 512])
ylim([0 512])

% Convert the format: Store them into individual nuclei traces
spot_fluo_traces = zeros(length(sch_pass),final_frame);
x_pos = zeros(length(sch_pass),final_frame);
y_pos = zeros(length(sch_pass),final_frame);

for i = 1:sch_final_num
    for j = 1:final_frame
        spot_fluo_traces(i,j) = processed_data(j).SpotFluoSum(i);
        x_pos(i,j) = processed_data(j).xcoord(i);
        y_pos(i,j) = processed_data(j).ycoord(i);
    end
end
%% Part 4: Process traces

% We would like to combine everything together and interpolate values

% Define time point to interpolate
time_point_final = -20:1/6:0;
x = ElapsedTime(1:final_frame)-ElapsedTime(final_frame);

spot_fluo_ins= zeros(sch_final_num,length(time_point_final));


for i = 1:sch_final_num
    
    spot_temp = spot_fluo_traces(i,:);

    % Interpolate traces
    spot_interp = interp1(x,spot_temp,time_point_final,'pchip');

    % pmoving average
    %spot_interp_avr = movmean(spot_interp,3);
    spot_interp_avr = movmean(spot_interp,1);
    spot_fluo_ins(i,:) = spot_interp_avr;
end

spot_fluo_full = [spot_fluo_full spot_fluo_ins'];

end

flag = sum(spot_fluo_full,1) > 0;
spot_fluo_full = spot_fluo_full(:,flag);

%% Plot combined late element dynamics

spot_fluo_mean = mean(spot_fluo_full,2);
spot_fluo_movmean = movmean(spot_fluo_mean,15);
spot_fluo_ste = std(spot_fluo_full,1,2)/sqrt(size(spot_fluo_full,2));
spot_fluo_ste_movmean = movmean(spot_fluo_ste,15);

fig = figure;
hold on
errorbar(time_point_final, spot_fluo_movmean, spot_fluo_ste_movmean);
%plot(time_point_final, spot_fluo_movmean);
%plot(time_point_final, spot_fluo_movmean/max(spot_fluo_movmean),'LineWidth',2);
xlim([-20 -2])
ylim([0 1E5])
%ylim([0 1.1])
pbaspect([3 2 1])

save('./data/late_dynamics_yellow.mat','time_point_final', 'spot_fluo_movmean', 'spot_fluo_ste_movmean')

