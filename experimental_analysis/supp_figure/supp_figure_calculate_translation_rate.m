clear all
close all
clc

%% Part 1: Read and display data

% embryo 1
load('./data/translation_rate/translation_rate_embryo1/translation_rate_embryo1_lin.mat');
load('./data/translation_rate/translation_rate_embryo1/CompiledParticles.mat')
final_frame = 82; % embryo 1
time_point_final = 0:1/6:29; % embryo 1


% embryo 2
%load('./data/translation_rate/translation_rate_embryo2/translation_rate_embryo2_lin.mat');
%load('./data/translation_rate/translation_rate_embryo2/CompiledParticles.mat')
%final_frame = 73; % embryo 2
%time_point_final = 0:1/6:23; % embryo 2


offset = 508892.34;
gamma_p = 1/7.9; % 1/min, protein half-life


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

% Convert the format: Store them into individual nuclei traces
nuclei_fluo_traces = zeros(length(sch_pass),final_frame);
sum_spot_fluo_traces = zeros(length(sch_pass),final_frame);
int_spot_fluo_traces = zeros(length(sch_pass),final_frame);
spot_fluo_traces = zeros(length(sch_pass),final_frame);

for i = 1:length(sch_pass)
    for j = 1:final_frame
        nuclei_fluo_traces(i,j) = processed_data(j).NuclearFluo(i);
        sum_spot_fluo_traces(i,j) = processed_data(j).SpotFluoSum(i);
        int_spot_fluo_traces(i,j) = processed_data(j).SpotFluoInt(i);
        spot_fluo_traces(i,j) = processed_data(j).SpotFluo(i);
    end
end


%% Part 4: Process traces
% We would like to combine everything together and interpolate values

% Define time point to interpolate
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
    nuclei_interp_avr = movmean(nuclei_interp,15);

    spot_fluo_ins(i,:) = spot_interp_avr;
    nuclei_fluo_ins(i,:) = nuclei_interp_avr;

    %spot_fluo_input(i,:) = movmean(spot_interp_avr,15);
    % use integrated mRNA level as input
    spot_fluo_input(i,:) = int_spot_interp;
    nuclei_fluo_output(i,:) = nuclei_interp_avr;

    spot_temp = spot_fluo_input(i,:);

end


%% Part 5: Estimate translation related rate (based on optimization)

ab_range = linspace(0,1,301); % combined rate for protein production
dt = 1/6; % min

ab_est = zeros(length(sch_pass),1);
ab_est_flag = zeros(length(sch_pass),1);

start_frame =  1;


for i = 1:length(sch_pass)
    int_spot_temp = spot_fluo_input(i,:);
    nuclei_temp = nuclei_fluo_output(i,:)-offset;
    
    int_protein_best = zeros(length(time_point_final),1);
    MSD_min = Inf;
    
    if sum(int_spot_temp) > 0
    
        % initialize array to store calculated MSD
        MSD = zeros(length(ab_range),1);
        
        for j = 1:length(ab_range)
            int_protein = zeros(length(time_point_final),1);
            int_protein(start_frame) = nuclei_temp(start_frame);

            for k = start_frame+1:length(time_point_final)
                int_protein(k) = int_protein(k-1)*2^(-dt*gamma_p) + ab_range(j)*int_spot_temp(k-1)*dt;
            end

            MSD(j) = sum((int_protein(start_frame:length(time_point_final))' - nuclei_temp(start_frame:length(time_point_final))).^2);
            if MSD(j)<MSD_min
                MSD_min = MSD(j);
                int_protein_best = int_protein;
            end
        end
        
        [temp, ind] = min(MSD);
        ab_est(i) = ab_range(ind);

        if max(int_spot_temp)>0.75E6

            ab_est_flag(i) = 1;

        end
    end
end

ab_final = ab_est((ab_est>0) & (ab_est<0.2) & (ab_est_flag==1));

%% bootstrap estimate
nBoots = 100;

boot_samples = bootstrp(nBoots,@nanmean,ab_final);
ab_final_mean = nanmean(boot_samples)
ab_final_ste = std(boot_samples)

%% plot figure

fig = figure;
histogram(ab_final,19,'Normalization','Probability')
title('Estimating \alpha')
xlabel('Estimated \alpha')
ylabel('Frequency')
xlim([0 0.2])
pbaspect([3 2 1])