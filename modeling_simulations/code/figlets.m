clc; clear; close all;

color_ftz_green = [122,167,116]/255;
color_autoreg_red = [212,108,85]/255;
color_early_blue = [114,141,192]/255;
color_latestart_gray = [112,112,112]/255;
color_shutdown_gold = [237,177,32]/255;

time_gastrulation = 20;

fi_bi = ftz_instance();
fi_mon_low = ftz_instance();
fi_mon_high = ftz_instance();

fi_mon_low.set('c',0.3);
fi_mon_high.set('n',0.7,'c',0.65);

PP = linspace(0,3e6,100);
tau_set = inf;

uic = [0,0];
a0_hi = 2e5;
a0_low = 1.25e5;
r0 = 1e4;
p0 = 2e5;


%% BISTABLE INTERSECTION PLOT (EMPIRICAL)
fi_bi.plot_ss_graphical_test(PP,tau_set,'protein','new_fig')
legend off;
title(''); xlabel('');
axis square;
xlim([0,max(PP)]);
ylim([0,fi_bi.gammaP*max(PP)]);
set(gcf,'position',[400 750 150 150]);


%% EXAMPLE BISTABLE INTERSECTION PLOTS
figure
subplot(1,3,1)
fi_mon_low.plot_ss_graphical_test(PP,tau_set)
legend off;
title(''); xlabel('');
axis square;
xlim([0,max(PP)]);
ylim([0,fi_bi.gammaP*max(PP)]);


subplot(1,3,2)
fi_bi.plot_ss_graphical_test(PP,tau_set)
legend off;
title(''); xlabel('');
axis square;
xlim([0,max(PP)]);
ylim([0,fi_bi.gammaP*max(PP)]);

subplot(1,3,3)
fi_mon_high.plot_ss_graphical_test(PP,tau_set)
legend off;
title(''); xlabel('');
axis square;
xlim([0,max(PP)]);
ylim([0,fi_bi.gammaP*max(PP)]);

set(gcf,'position',[550 750 550 180])


%% VISUALIZE EARLY ELEMENT BEHAVIOR
paxlim = 1.2e6;
scc = 5;    % scale constant for visualizing decay better

tarr = linspace(0,140,100);
fi_bi.set('a0',a0_hi,'r0',3e5,'p0',2e4);
[~,P_bi_high] = fi_bi.simulate(uic,tarr,'no_plot','no_save','off');

figure

% lines for "half-life" of production rate decay
plot((1/fi_bi.b)*[1,1]-time_gastrulation,scc*fi_bi.c1*exp(-1)*[0,1],'--','color',0.7*[1 1 1]); hold on;
plot((1/fi_bi.b)*[0,1]-time_gastrulation,scc*fi_bi.c1*exp(-1)*[1,1],'--','color',0.7*[1 1 1]);

hs = plot(tarr-time_gastrulation,scc*fi_bi.c1*exp(-fi_bi.b*tarr), ...
    tarr-time_gastrulation,fi_bi.r(tarr), ...
    tarr-time_gastrulation,fi_bi.p(tarr));
hold off;
xlim([-time_gastrulation tarr(end)-time_gastrulation])
ylim([0,paxlim]); yticks([]);
xlabel('time (min)')
set_figure_defaults(gcf)
set(gcf,'position',[560 768 550 180])

legend(hs,'r(t)','R(t)','P(t)')


%% VISUALIZE EARLY ELEMENT SHUTOFF
fi_delay = ftz_instance();
fi_delay.set('a0',a0_hi,'r0',r0,'p0',p0);

paxlim = 1.2e6;

sd_time = 10;

tarr = linspace(0,40,200);
[~,ix_to_sd] = min(abs(tarr - sd_time));

a1 = nan(numel(tarr),1);
r1 = nan(numel(tarr),1);
p1 = nan(numel(tarr),1);

a1(1:ix_to_sd) = fi_delay.c1*exp(-fi_delay.b*tarr(1:ix_to_sd));
r1(1:ix_to_sd) = fi_delay.r(tarr(1:ix_to_sd));
p1(1:ix_to_sd) = fi_delay.p(tarr(1:ix_to_sd));

fi_delay.set('a0',0,'r0',r1(ix_to_sd),'p0',p1(ix_to_sd));
a1(ix_to_sd+1:end) = 0;
r1(ix_to_sd+1:end) = fi_delay.r(tarr(ix_to_sd+1:end) - tarr(ix_to_sd));
p1(ix_to_sd+1:end) = fi_delay.p(tarr(ix_to_sd+1:end) - tarr(ix_to_sd));

figure
plot((sd_time-time_gastrulation)*[1 1],[0,paxlim],'--','color','b');
hold on;
plot([0,0],[0,paxlim],'--','color',0.7*[1 1 1]);
hs = plot(tarr - time_gastrulation,a1,tarr - time_gastrulation,r1,tarr - time_gastrulation,p1);
hold off;

legend(hs,'a(t)','r(t)','p(t)')
ylim([0,paxlim])
set_figure_defaults(gcf)
set(gcf,'position',[560 768 550 180])


%% TOFF VS TON (TOFF FIRST)
paxlim = 1.5e6;
tarr = linspace(0,140,1000);

sd_time = 10;
ls_time = 15;
[~,ix_to_sd] = min(abs(tarr - sd_time));
[~,ix_to_ls] = min(abs(tarr - ls_time));
fi_delay.set('a0',a0_hi,'r0',1e6,'p0',1e6);

adelay = zeros(numel(tarr),1);
rdelay = nan(numel(tarr),1);
pdelay = nan(numel(tarr),1);
Rdelay = zeros(numel(tarr),1);
Pdelay = zeros(numel(tarr),1);

adelay(1:ix_to_sd) = fi_delay.c1*exp(-fi_delay.b*tarr(1:ix_to_sd));
rdelay(1:ix_to_sd) = fi_delay.r(tarr(1:ix_to_sd));
pdelay(1:ix_to_sd) = fi_delay.p(tarr(1:ix_to_sd));

fi_delay.set('a0',0,'r0',rdelay(ix_to_sd),'p0',pdelay(ix_to_sd));

rdelay(ix_to_sd:ix_to_ls) = fi_delay.r(tarr(ix_to_sd:ix_to_ls) - tarr(ix_to_sd));
pdelay(ix_to_sd:ix_to_ls) = fi_delay.p(tarr(ix_to_sd:ix_to_ls) - tarr(ix_to_sd));

fi_delay.set('r0',rdelay(ix_to_ls),'p0',pdelay(ix_to_ls));

[R,P] = fi_delay.simulate([0,0],tarr(ix_to_ls:end) - tarr(ix_to_ls),'no_plot','no_save');
rdelay(ix_to_ls:end) = fi_delay.r(tarr(ix_to_ls:end) - tarr(ix_to_ls))';
pdelay(ix_to_ls:end) = fi_delay.p(tarr(ix_to_ls:end) - tarr(ix_to_ls))';
Rdelay(ix_to_ls:end) = R;
Pdelay(ix_to_ls:end) = P;

figure
subplot(4,1,1)
plot(tarr - time_gastrulation,adelay,'color',color_shutdown_gold);
xlim([tarr(1),tarr(end)] - time_gastrulation);
xticks([])
yticks([])

subplot(4,1,2)
plot(tarr - time_gastrulation,(tarr > ls_time),'color',color_latestart_gray);
xlim([tarr(1),tarr(end)] - time_gastrulation);
xticks([])
yticks([])

subplot(4,1,3:4)
plot((sd_time-time_gastrulation)*[1 1],[0,paxlim],'--','color',color_early_blue);
hold on;
plot((ls_time-time_gastrulation)*[1 1],[0,paxlim],'--','color',color_autoreg_red);
hs = plot(tarr - time_gastrulation, pdelay, 'color', color_early_blue);
hs = [hs,plot(tarr - time_gastrulation, Pdelay, 'color', color_autoreg_red)];
hs = [hs,plot(tarr - time_gastrulation, pdelay + Pdelay, 'color', color_ftz_green)];
hold off;

legend(hs,'early Ftz','late Ftz','total Ftz')
xlim([tarr(1),tarr(end)] - time_gastrulation);
ylim([0,paxlim]);
yticks([]);
% xlabel('time (min)')
set_figure_defaults(gcf);
set(gcf,'position',[560 768 375 200])


%% TOFF VS. TON (TON FIRST)
paxlim = 2e6;
tarr = linspace(0,140,1000);

sd_time = 15;
ls_time = 10;
[~,ix_to_sd] = min(abs(tarr - sd_time));
[~,ix_to_ls] = min(abs(tarr - ls_time));
fi_delay.set('a0',a0_hi,'r0',0.8e6,'p0',1e6);

adelay = zeros(numel(tarr),1);
rdelay = nan(numel(tarr),1);
pdelay = nan(numel(tarr),1);
Rdelay = zeros(numel(tarr),1);
Pdelay = zeros(numel(tarr),1);

adelay(1:ix_to_sd) = fi_delay.c1*exp(-fi_delay.b*tarr(1:ix_to_sd));
rdelay(1:ix_to_sd) = fi_delay.r(tarr(1:ix_to_sd));
pdelay(1:ix_to_sd) = fi_delay.p(tarr(1:ix_to_sd));

fi_delay.set('a0',adelay(ix_to_ls),'r0',rdelay(ix_to_ls), ...
    'p0',pdelay(ix_to_ls));
[R,P] = fi_delay.simulate([0,0],tarr(ix_to_ls:ix_to_sd) - tarr(ix_to_ls), ...
    'no_plot','no_save');
Rdelay(ix_to_ls:ix_to_sd) = R;
Pdelay(ix_to_ls:ix_to_sd) = P;

fi_delay.set('a0',0,'r0',rdelay(ix_to_sd),'p0',pdelay(ix_to_sd));
[R,P] = fi_delay.simulate([Rdelay(ix_to_sd),Pdelay(ix_to_sd)], ...
    tarr(ix_to_sd:end) - tarr(ix_to_sd),'no_plot','no_save');
rdelay(ix_to_sd:end) = fi_delay.r(tarr(ix_to_sd:end)-tarr(ix_to_sd));
pdelay(ix_to_sd:end) = fi_delay.p(tarr(ix_to_sd:end)-tarr(ix_to_sd));
Rdelay(ix_to_sd:end) = R;
Pdelay(ix_to_sd:end) = P;


figure
subplot(4,1,1)
plot(tarr - time_gastrulation,(tarr > ls_time),'color',color_latestart_gray);
xlim([tarr(1),tarr(end)] - time_gastrulation);
xticks([])
yticks([])

subplot(4,1,2)
plot(tarr - time_gastrulation,adelay,'color',color_shutdown_gold);
xlim([tarr(1),tarr(end)] - time_gastrulation);
xticks([])
yticks([])

subplot(4,1,3:4)
plot((sd_time-time_gastrulation)*[1 1],[0,paxlim],'--','color',color_shutdown_gold);
hold on;
plot((ls_time-time_gastrulation)*[1 1],[0,paxlim],'--','color',color_latestart_gray);
% plot([0,0],[0,paxlim],'--','color',0.7*[1 1 1]);
% hs = plot(tarr - time_gastrulation,Rdelay,tarr - time_gastrulation,Pdelay);
hs = plot(tarr - time_gastrulation, pdelay, 'color', color_early_blue);
hs = [hs,plot(tarr - time_gastrulation, Pdelay, 'color', color_autoreg_red)];
hs = [hs,plot(tarr - time_gastrulation, pdelay + Pdelay, 'color', color_ftz_green)];
hold off;

% legend(hs,'R_{tot}(t)','P_{tot}(t)')
legend(hs,'early Ftz','late Ftz','total Ftz')
xlim([tarr(1),tarr(end)] - time_gastrulation);
ylim([0,paxlim]);
yticks([]);
% xlabel('time (min)')
set_figure_defaults(gcf);
set(gcf,'position',[560 768 375 200])