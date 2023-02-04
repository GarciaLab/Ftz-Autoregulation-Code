clear all
close all
clc

addpath(genpath('../lib'))

%% Read data
early_dynamics = load('data/early_dynamics.mat');
late_dynamics = load('data/late_dynamics.mat');

%% Plot data

%norm_factor_early = max(early_dynamics.spot_fluo_movmean);
%norm_factor_late = max(late_dynamics.spot_fluo_movmean);

norm_factor_early = early_dynamics.spot_fluo_movmean(1); % 1, -20min
norm_factor_late = late_dynamics.spot_fluo_movmean(109); % 109, -2min

fig = figure;
hold on
yyaxis left
plot(early_dynamics.time_point_final,early_dynamics.spot_fluo_movmean/norm_factor_early,'LineWidth',2);
boundedline(early_dynamics.time_point_final,early_dynamics.spot_fluo_movmean/norm_factor_early,early_dynamics.spot_fluo_ste_movmean/norm_factor_early,'-','nan', 'gap','alpha');
ylabel('early element level (au)')
ylim([0 1.25])

yyaxis right
plot(late_dynamics.time_point_final,late_dynamics.spot_fluo_movmean/norm_factor_late,'LineWidth',2);
boundedline(late_dynamics.time_point_final,late_dynamics.spot_fluo_movmean/norm_factor_late,late_dynamics.spot_fluo_ste_movmean/norm_factor_late,'-','nan', 'gap','alpha');
ylabel('late element level (au)')
ylim([0 1.25])

xlim([-20 -2])
xlabel('time (min)')

pbaspect([3 2 1])


%% Plot original data (not working)

%correction_factor = 1.75; % 
correction_factor = 15/10*0.87/0.64;

fig = figure;
hold on
plot(early_dynamics.time_point_final,early_dynamics.spot_fluo_movmean,'LineWidth',2);
boundedline(early_dynamics.time_point_final,early_dynamics.spot_fluo_movmean,early_dynamics.spot_fluo_ste_movmean,'-','nan', 'gap','alpha');
ylabel('early element level (au)')


plot(late_dynamics.time_point_final,late_dynamics.spot_fluo_movmean*correction_factor,'LineWidth',2);
boundedline(late_dynamics.time_point_final,late_dynamics.spot_fluo_movmean*correction_factor,late_dynamics.spot_fluo_ste_movmean*correction_factor,'-','nan', 'gap','alpha');
ylabel('late element level (au)')

xlim([-20 -2])
xlabel('time (min)')
ylim([0 1.25E5])

pbaspect([3 2 1])
