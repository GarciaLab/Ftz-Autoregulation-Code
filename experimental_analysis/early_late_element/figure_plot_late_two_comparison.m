clear all
close all
clc

addpath(genpath('../lib'))

%% Read data
late_dynamics = load('data/late_dynamics.mat');
late_dynamics_yellow = load('data/late_dynamics_yellow.mat');

%% Plot data

norm_factor_late = late_dynamics.spot_fluo_movmean(109); % 109, -2min
norm_factor_late_yellow = late_dynamics_yellow.spot_fluo_movmean(109); % 109, -2min

fig = figure;
hold on
yyaxis left
plot(late_dynamics.time_point_final,late_dynamics.spot_fluo_movmean/norm_factor_late,'LineWidth',2);
boundedline(late_dynamics.time_point_final,late_dynamics.spot_fluo_movmean/norm_factor_late,late_dynamics.spot_fluo_ste_movmean/norm_factor_late,'-','nan', 'gap','alpha');
ylabel('late element level (au)')
ylim([0 1.25])

yyaxis right
plot(late_dynamics_yellow.time_point_final,late_dynamics_yellow.spot_fluo_movmean/norm_factor_late_yellow,'LineWidth',2);
boundedline(late_dynamics_yellow.time_point_final,late_dynamics_yellow.spot_fluo_movmean/norm_factor_late_yellow,late_dynamics_yellow.spot_fluo_ste_movmean/norm_factor_late_yellow,'-','nan', 'gap','alpha');
ylabel('late element (yellow) level (au)')
ylim([0 1.25])

xlim([-20 -2])
xlabel('time (min)')

pbaspect([3 2 1])