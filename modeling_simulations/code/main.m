clc; clear; close all;

% -- CHOOSE WHICH PARTS OF SCRIPT TO RUN -- %
analysis_option = 'prediction_lite';
% options are:
% - 'estimate_c'
% - 'prediction'
% - 'prediction_lite' (binary classification only)
% - 'stochastic_parameter_estimation' (mRNA/protein noise)
% - 'stochastic_simulation'

% --- %

color_ftz_green = [122,167,116]/255;
color_autoreg_red = [212,108,85]/255;
color_early_blue = [114,141,192]/255;

color_green_hc = [153 255 0]/255;
color_latestart_gray = [112,112,112]/255;
color_shutdown_gold = [237,177,32]/255;

ssf = switching_separatrix('switching_separatrix_a0-r0-p0.mat');
fi = ssf.fi;
tssf = switching_separatrix('trans_switching_separatrix_a0-r0-p0.mat');

embryoid_to_ignore = [];
rot_emb = [];   % -1 to rotate 180 deg


% parts of the script to run
estimate_c = false;
simulate_trajectories = false;
estimate_gap = false;   % requires simulate_trajectories to be true
estimate_window = false; % requires simulate_trajectories to be true
get_noise_params = false;
save_noise_params = false;
stochastic_sims = false;
stripe_to_analyze = 4;
embryo_set = 'prediction';

switch analysis_option
    case 'estimate_c'
        embryo_set = 'measurement';
        estimate_c = true;
    case 'prediction'
        simulate_trajectories = true;
        estimate_gap = true;
        estimate_window = true;
    case 'prediction_lite'
        simulate_trajectories = true;
    case 'stochastic_parameter_estimation'
        embryo_set = 'stochastic_param';
        get_noise_params = true;
        % save_noise_params = true;
    case 'stochastic_simulation'
        embryo_set = 'stochastic';
        stochastic_sims = true;
end

switch embryo_set
    case 'measurement'
        load('pooled_embryos_io_boundary_smooth_traces.mat');
        rot_emb = [1 -1 1 1 1 1 -1];
    case 'prediction'
        % ignore embryoid 1 (embryo 6) which appears to be in a weird D/V orientation
        load('pooled_embryos_excluded_io_boundary_smooth_traces.mat');
        rot_emb = [-1,-1,1,-1]; embryoid_to_ignore = 1;
    case 'stochastic_param'
        load('pooled_embryos_io_all_smooth_traces.mat');
        rot_emb = [1 -1 1 1 1 1 -1];
        stripe_to_analyze = 'all';
    case 'stochastic'
        load('pooled_embryos_excluded_io_boundary_smooth_traces.mat');
        rot_emb = [-1,-1,1,-1];
        embryoid_to_ignore = [1,3]; % 1
end


nix_per_min = 6;    % number of indices per minute
nix_for_avg = 6;
time_autoreg_start = -20;	% min
time_sim_start = -20;   % min

[~,ix_autoreg_start] = min(abs(time_point_final - time_autoreg_start));
[~,ix_gastrulation] = min(abs(time_point_final));
[~,ix_sim_start] = min(abs(time_point_final - time_sim_start));

begin_ix_for_avg = (ix_autoreg_start-nix_for_avg):ix_autoreg_start;
end_ix_for_avg = (time_ix(end)-nix_for_avg:time_ix(end));

% on/off cutoff
prot_on_thresh = fi.prot_on_thresh;

pred_by_sim = true;
stochastic_sim = false;
plot_select_nuc = true;

time_dependent_iof = false;



%% ESTIMATE EARLY VS. LATE RNA AND PROTEIN
valid_indices_set = false;
prot = fillmissing(prot,'linear');


a0 = nanmean(((mRNA(begin_ix_for_avg+1,:) - ...
    mRNA(begin_ix_for_avg,:)) / dt) + ...
    fi.gammaR*mRNA(begin_ix_for_avg,:),1);   % last term corrects for decay
a0(a0 < 0) = 0;

r0 = nanmean(mRNA(begin_ix_for_avg,:),1);
p0 = nanmean(prot(begin_ix_for_avg,:),1);

mRNAearly_est = nan(numel(time_point_final),nnuclei);
Pearly_est = nan(numel(time_point_final),nnuclei);
a = nan(numel(time_point_final),nnuclei);
for ii = 1:nnuclei
    fi.set('a0',a0(ii),'r0',r0(ii),'p0',p0(ii));
    mRNAearly_est(:,ii) = fi.r(time_point_final-time_autoreg_start);
    Pearly_est(:,ii) = fi.p(time_point_final-time_autoreg_start);
    a(:,ii) = fi.c1*exp(-fi.b*(time_point_final-time_autoreg_start))';
end

Plate_est = 0*Pearly_est;
Plate_est(ix_autoreg_start+1:end,:) = prot(ix_autoreg_start+1:end,:) - ...
    Pearly_est(ix_autoreg_start+1:end,:);
Plate_est(Plate_est < 0) = 0;

mRNAlate_est = 0*mRNAearly_est;
mRNAlate_est(ix_autoreg_start+1:end,:) = ...
    fillmissing(mRNA(ix_autoreg_start+1:end,:),'linear') - mRNAearly_est(ix_autoreg_start+1:end,:);
mRNAlate_est(mRNAlate_est < 0) = 0;


%% CALCULATE a0 AND SET VALID INDICES
if ~valid_indices_set
    mRNA_final = nanmean(mRNA(end_ix_for_avg,:),1);
    prot_final = nanmean(prot(end_ix_for_avg,:),1);
    
    valid_ix = ~ismember(embryoid',embryoid_to_ignore);
    for cur_embryoid = 1:nembryos
        cur_emb_ix = (embryoid == cur_embryoid);
        x_midpoint = mean([min(rot_emb(cur_embryoid)*xcoord(cur_emb_ix)), ...
            max(rot_emb(cur_embryoid)*xcoord(cur_emb_ix))]);
        switch stripe_to_analyze
            case 'all'
            case 3
                valid_ix(cur_emb_ix) = ...
                    valid_ix(cur_emb_ix) & (rot_emb(cur_embryoid)*xcoord(cur_emb_ix)' < x_midpoint);
            case 4
                valid_ix(cur_emb_ix) = ...
                    valid_ix(cur_emb_ix) & (rot_emb(cur_embryoid)*xcoord(cur_emb_ix)' > x_midpoint);
            otherwise
                error(['no data to analyze stripe ' num2str(stripe_to_analyze)])
        end
    end
    
    nvalid_nuc = sum(valid_ix);
    
    a0 = a0(valid_ix);
    mRNA = mRNA(:,valid_ix);
    mRNA_final = mRNA_final(valid_ix);
    prot = prot(:,valid_ix);
    prot_final = prot_final(valid_ix);
    protein_trace = protein_trace(valid_ix,:);
    
    a = a(:,valid_ix);
    Plate = Plate_est(:,valid_ix);
    Pearly = Pearly_est(:,valid_ix);
    
    mRNAlate = mRNAlate_est(:,valid_ix);
    mRNAearly = mRNAearly_est(:,valid_ix);
    
    valid_indices_set = true;
end

ix_on = (prot_final > prot_on_thresh);


%% ESTIMATE C FROM ANALYTICAL SOLUTION
% estimate c from analytical solution to a linear system model for protein
% level at gastrulation

% for time_cest_start = -3, c_est_final = 0.4491
if estimate_c
    Pfinthresh = 1.5e6;
    
    % personal observations of accuracy from forward simulation suggest
    % the time window should be no larger than 2.5 min
    %
    % choose closer to gastrulation because gene regulatory function
    % expected to be closest to final shape
    time_cur_cest_start = -3;	% min
    time_cur_cest_end = 0;      % fix at gastrulation b/c f is highest
    [~,ix_cest_start] = min(abs(time_point_final - time_cur_cest_start));
    [~,ix_cest_end] = min(abs(time_point_final - time_cur_cest_end));
    
    if ix_cest_end == ix_gastrulation
        % correction because we only have an mRNA estimate up until the
        % time point before gastrulation
        ix_cest_end = ix_cest_end - 1;
        time_cur_cest_end = time_cur_cest_end - (1/nix_per_min);
    end
    
    % pick nuclei with close to maximum activity from early element from
    % start of autoregulation to gastrulation (necessary for assuming
    % linear dynamics)
    ixs2use = (prot(ix_cest_end,:) > Pfinthresh) & (prot(ix_cest_start,:) > Pfinthresh);
    
    tfin = time_cur_cest_end - time_cur_cest_start;
    c_all_traces = nan(numel(ixs2use),1);
    
    for ii = 1:numel(ixs2use)
        fi_new = ftz_instance();
        if ixs2use(ii)
            fi_new.set('a0',a(ix_cest_start,ii), ...
                'r0',mRNAearly(ix_cest_start,ii), ...
                'p0',Pearly(ix_cest_start,ii));
            
            c_all_traces(ii) = ...
                (fi_new.gammaP/(fi_new.f(inf,Pfinthresh)*(1-exp(-fi_new.gammaP*tfin)))) * ...
                ((mRNAlate(ix_cest_end,ii) - mRNAlate(ix_cest_start,ii)*exp(-fi.gammaP*tfin)) - ...
                ((fi_new.gammaP-fi_new.gammaR)/fi_new.alpha) * ...
                (Plate(ix_cest_end,ii) - Plate(ix_cest_start,ii)*exp(-fi.gammaP*tfin)));
        end
    end
    c_est = nanmean(c_all_traces);
    c_est_stderr = nanstd(c_all_traces)/sqrt(sum(~isnan(c_all_traces)));
    disp(['estimated c from analytical solution = ' num2str(c_est)]);
end


%% ESTIMATE C FROM ANALYTICAL SOLUTION (SLIDING WINDOW)
% estimate c from analytical solution to a linear system model for protein
% level at gastrulation

% for time_cest_start = -3, c_est_final = 0.4454
if estimate_c
    Pfinthresh = 1.5e6;
    
    % average over estimates for sliding window
    nmin_per_window = 1/nix_per_min;
    
    % short time window before gastrulation
    time_cest_start = -3;
    time_cest_end = 0;
    
    [~,ix_cest_start] = min(abs(time_point_final - time_cest_start));
    [~,ix_cest_end] = min(abs(time_point_final - time_cest_end));
    
    if ix_cest_end == ix_gastrulation
        % correction because we only have an mRNA estimate up until the
        % time point before gastrulation
        ix_cest_end = ix_cest_end - 1;
        time_cest_end = time_cest_end - (1/nix_per_min);
    end
    
    % pick nuclei with close to maximum activity from early element from
    % start of autoregulation to gastrulation (necessary for assuming
    % linear dynamics)
    ixs2use = (prot(ix_cest_end,:) > Pfinthresh) & (prot(ix_cest_start,:) > Pfinthresh);
    
    ix_per_win = nmin_per_window*nix_per_min;
    tfin = nmin_per_window;
    
    nwin = ix_cest_end - ix_cest_start - ix_per_win + 1;
    c_est = nan(nwin,1);
    c_err = nan(nwin,1);
    win_start_ixs = ix_cest_start:(ix_cest_end-ix_per_win);
    c_all_traces = nan(nwin,numel(ixs2use));
    for ii_win = 1:nwin
        ix_win = win_start_ixs(ii_win);
        ix_win_end = ix_win + ix_per_win;
        
        for ii = 1:numel(ixs2use)
            fi_new = ftz_instance();
            if ixs2use(ii)
                fi_new.set('a0',a(ix_win,ii), ...
                    'r0',mRNAearly(ix_win,ii), ...
                    'p0',Pearly(ix_win,ii));
                
                c_all_traces(ii_win,ii) = ...
                    (fi_new.gammaP/(fi_new.f(inf,Pfinthresh)*(1-exp(-fi_new.gammaP*tfin)))) * ...
                    ((mRNAlate(ix_win_end,ii) - mRNAlate(ix_win,ii)*exp(-fi.gammaP*tfin)) - ...
                    ((fi_new.gammaP-fi_new.gammaR)/fi_new.alpha) * ...
                    (Plate(ix_win_end,ii) - Plate(ix_win,ii)*exp(-fi.gammaP*tfin)));
            end
        end
    end
    
    c_est_final = nanmean(c_all_traces(:));
    c_stderr_final = nanstd(c_all_traces(:))/sqrt(sum(~isnan(c_all_traces(:))));
    disp(['estimated c from analytical solution with sliding window = ' ...
        num2str(c_est_final) ' +/- ' num2str(c_stderr_final)]);
    
    
    %% SUPPLEMENT FIGURES
    Prange = linspace(0,3e6,100);
    prot2plot = prot(ix_cest_start:ix_cest_end,ixs2use);
    prot2plot = prot2plot(:);
    sf_temp = spot_fluo_output(valid_ix,:);
    autoreg2plot = sf_temp(ixs2use,ix_cest_start:ix_cest_end)';
    autoreg2plot = autoreg2plot(:);
    
    protbelowthresh = prot(ix_cest_start:ix_cest_end,~ixs2use);
    protbelowthresh = protbelowthresh(:);
    autoregbelowthresh = sf_temp(~ixs2use,ix_cest_start:ix_cest_end)';
    autoregbelowthresh = autoregbelowthresh(:);
    
    figure
    scatter(prot2plot,autoreg2plot,[],0.8*[1,1,1],'filled')
    hold on
    scatter(protbelowthresh,autoregbelowthresh,[],0.9*[1,1,1],'filled')
    plot(Prange,fi.f(0,Prange),'color',color_autoreg_red)
    plot(Pfinthresh*[1,1],[0,1.2e6],'k--')
    hold off
    ylim([0,1.2e6]);
    xlabel('[Ftz] (au)')
    ylabel('mRNA production rate')
    set_figure_defaults(gcf)
    
    
    %find(ixs2use) % choose indices from here
    sample_trace_ix = [28,112,182,249];
    
    figure
    subplot(3,1,1)
    plot(time_point_final,mRNAlate(:,sample_trace_ix));
    xlim([-10,0])
    xlabel('time (min)')
    ylabel('R_{late}(t)')
    
    tprime = -3;
    tdelt = nmin_per_window;
    subplot(3,1,2)
    plot(time_point_final,Plate(:,sample_trace_ix));
    hold on;
    plot(tprime*[1,1],[0,1.2e6],'k--')
    plot((tprime+tdelt)*[1,1],[0,1.2e6],'k--')
    hold off;
    xlim([-10,0])
    xlabel('time (min)')
    ylabel('P_{late}(t)')
    ylim([0,1e6]);
    
    subplot(3,1,3)
    plot(time_point_final(win_start_ixs),c_all_traces,'color',[0 0 0 0.1])
    hold on;
    plot(time_point_final(win_start_ixs),nanmean(c_all_traces,2),'color',color_autoreg_red)
    hold off;
    xlabel('\tau (min)')
    ylabel('c')
    
    set_figure_defaults(gcf)
    set(gcf,'position',[550 325 550 650])
    
    figure
    histogram(c_all_traces(~isnan(c_all_traces)))
    xlabel('c')
    
end


%% SIMULATE TRAJECTORIES
if simulate_trajectories
    if time_dependent_iof
        cur_a = fi.a;
        cur_n = fi.n;
        cur_K = fi.K;
        
        tau_half = 5;
        beta = log(2)/tau_half;  % growth rate assuming iof reaches 50% of height in time tau_half
        t_on = 0;                % minutes past simulation start
        
        iof = @(t,x) (t > t_on)*(1-exp(-beta*(t-t_on)))* ...
            ((cur_a*x.^cur_n)./(cur_K.^cur_n + x.^cur_n));
        
        fi.set('f',iof);
    end
    
    tmax = 200-time_sim_start;
    tarr = 0:dt:tmax;
    
    aearly_all = nan(numel(tarr),nvalid_nuc);
    Rearly_all = nan(numel(tarr),nvalid_nuc);
    Rlate_all = nan(numel(tarr),nvalid_nuc);
    Pearly_all = nan(numel(tarr),nvalid_nuc);
    Plate_all = nan(numel(tarr),nvalid_nuc);
    Ptot_all = nan(numel(tarr),nvalid_nuc);
    for ii = 1:nvalid_nuc
        fi.set('a0',a(ix_sim_start,ii), ...
            'r0',mRNAearly(ix_sim_start,ii), ...
            'p0',Pearly(ix_sim_start,ii));
        
        uic = [mRNAlate(ix_sim_start,ii),Plate(ix_sim_start,ii)];
        [M,P] = fi.simulate(uic,tarr, ...
            'no_plot','no_save','off');
        
        aearly_all(:,ii) = fi.c1*exp(-fi.b*tarr)';
        Rearly_all(:,ii) = (1/(fi.gammaR-fi.b))*aearly_all(:,ii) + fi.c2*exp(-fi.gammaR*tarr)';
        Rlate_all(:,ii) = M;
        Pearly_all(:,ii) = fi.p(tarr)';
        Plate_all(:,ii) = P;
        Ptot_all(:,ii) = Pearly_all(:,ii) + P;
    end
    Rtot_all = Rearly_all + Rlate_all;
    
    time_point_sim = tarr + time_sim_start;
    [~,ix_mirrorautoreg] = min(abs(time_point_sim + time_autoreg_start));
    
    pred_on_thresh_sim = Ptot_all > prot_on_thresh;
    pred_off_thresh_sim = ~pred_on_thresh_sim;
    correct_pred_thresh_sim = (pred_on_thresh_sim & ix_on) | ...
        (pred_off_thresh_sim & ~ix_on);
    frac_correct_pred_thresh_sim = sum(correct_pred_thresh_sim/size(Ptot_all,2),2);
    
    pred_on_thresh_no_autoreg = Pearly_all > prot_on_thresh;
    pred_off_thresh_no_autoreg = ~pred_on_thresh_no_autoreg;
    correct_pred_thresh_no_autoreg = (pred_on_thresh_no_autoreg & ix_on) | ...
        (pred_off_thresh_no_autoreg & ~ix_on);
    frac_correct_pred_thresh_no_autoreg = sum(correct_pred_thresh_no_autoreg/size(Pearly_all,2),2);
    
    pred_on_thresh_emp = prot > prot_on_thresh;
    pred_off_thresh_emp = ~pred_on_thresh_emp;
    correct_pred_thresh_emp = (pred_on_thresh_emp & ix_on) | ...
        (pred_off_thresh_emp & ~ix_on);
    frac_correct_pred_thresh_emp = sum(correct_pred_thresh_emp/size(Ptot_all,2),2);
    
    disp(['From direct simulation: ' ...
        num2str(100*frac_correct_pred_thresh_sim(ix_gastrulation)) '%'])
    disp(['From direct simulation, no autoregulation: ' ...
        num2str(100*frac_correct_pred_thresh_no_autoreg(ix_gastrulation)) '%']);
      
    
    %% PREDICT FROM SWITCHING SEPARATRIX
    if time_dependent_iof
        warning(['time-dependent iof was chosen for simulations, but ' ...
            'switching separatrices are computed for time-independent iof']);
    end
    
    if ix_sim_start > ix_autoreg_start
        warning(['simulation begins after start of autoregulation, ' ...
            'so i.c.s for simulated trajectories are nonzero => ' ...
            'separatrix is a conservative bound'])
    end
    
    % use surface normal to determine if above or below surface
    % NOTE: surface has finite granularity so some data points very close
    % to the surface may not be correctly predicted due to numerical error
    fprintf('Steady state separatrix: ')
    ssf.prediction_accuracy(a(ix_sim_start,:),mRNAearly(ix_sim_start,:), ...
        Pearly(ix_sim_start,:),ix_on,valid_ix,embryoid);
    
    fprintf('Transient separatrix: ')
    tssf.prediction_accuracy(a(ix_sim_start,:),mRNAearly(ix_sim_start,:), ...
        Pearly(ix_sim_start,:),ix_on,valid_ix,embryoid);
    
    ix_disagree = xor(ssf.predicted_on,tssf.predicted_on);
    
    nuc_id = [27,16,84];
    nnuc_id = numel(nuc_id);
    
    a0_for_plot = a(ix_sim_start,:);
    mRNA_for_plot = mRNAearly(ix_sim_start,:);
    prot_for_plot = Pearly(ix_sim_start,:);
    
    figure
    h = ssf.surf();
    hold on;
    ht = tssf.surf();
    set(h,'facealpha',0.5,'facecolor',color_autoreg_red);
    set(ht,'facealpha',0.5,'facecolor',color_early_blue);
    xlabel('early mRNA production rate')
    ylabel('early mRNA (a.u.)')
    zlabel('early protein (a.u.)')
    
    % color by permanence of predicted fate
    h1 = scatter3(a0_for_plot(~tssf.predicted_on & ~tssf.false_negatives), ...
        mRNA_for_plot(~tssf.predicted_on & ~tssf.false_negatives), ...
        prot_for_plot(~tssf.predicted_on & ~tssf.false_negatives),75,'k','filled');
    h2 = scatter3(a0_for_plot(ix_disagree & ~tssf.false_positives), ...
        mRNA_for_plot(ix_disagree & ~tssf.false_positives), ...
        prot_for_plot(ix_disagree & ~tssf.false_positives),75, ...
        'k','filled');
    h3 = scatter3(a0_for_plot(ssf.predicted_on), ...
        mRNA_for_plot(ssf.predicted_on), ...
        prot_for_plot(ssf.predicted_on),75,0.65*color_green_hc,'filled');
    h4 = scatter3(a0_for_plot(tssf.false_positives), ...
        mRNA_for_plot(tssf.false_positives), ...
        prot_for_plot(tssf.false_positives),75, ...
        'r^','filled');
    h5 = scatter3(a0_for_plot(tssf.false_negatives), ...
        mRNA_for_plot(tssf.false_negatives), ...
        prot_for_plot(tssf.false_negatives),75, ...
        'rv','filled');
    
    view(34,5);
    
    hold off;
    
    l = legend([h1,h2,h3,h4,h5],{'fated off','transiently high','fated on', ...
        'false negative','false positive'},'location','northeast');
    cur_pos = get(gcf,'position');
    set(gcf,'position',[cur_pos(1),cur_pos(2),500,400]);
    set_figure_defaults(gcf,'font')
    set(l,'box','on')
    
    xlim([0,5e5])
    zlim([0,2.5e6])
    
    
    %% SEPARATE PLOT FOR TRANSIENT FATE
    nuc_id = [27,16,84];
    nnuc_id = numel(nuc_id);
    
    figure
    ht = tssf.surf(); hold on;
    set(ht,'facealpha',0.5,'facecolor',color_early_blue);
    xlabel('early mRNA production rate')
    ylabel('early mRNA (a.u.)')
    zlabel('early protein (a.u.)')
    
    pforcol = prot(ix_gastrulation,:)';
    pforcol = (pforcol - min(pforcol)) / range(pforcol);
    cols = pforcol*color_green_hc;
    
    % color by expression level at gastrulation
    h1 = scatter3(a0_for_plot(~tssf.predicted_on & ~tssf.false_negatives), ...
        mRNA_for_plot(~tssf.predicted_on & ~tssf.false_negatives), ...
        prot_for_plot(~tssf.predicted_on & ~tssf.false_negatives),75, ...
        cols(~tssf.predicted_on & ~tssf.false_negatives,:),'filled');
    h2 = scatter3(a0_for_plot(tssf.predicted_on & ~tssf.false_positives), ...
        mRNA_for_plot(tssf.predicted_on & ~tssf.false_positives), ...
        prot_for_plot(tssf.predicted_on & ~tssf.false_positives),75, ...
        cols(tssf.predicted_on & ~tssf.false_positives,:),'filled');
    h3 = scatter3(a0_for_plot(tssf.false_positives), ...
        mRNA_for_plot(tssf.false_positives), ...
        prot_for_plot(tssf.false_positives),75, ...
        'r^','filled');
    h4 = scatter3(a0_for_plot(tssf.false_negatives), ...
        mRNA_for_plot(tssf.false_negatives), ...
        prot_for_plot(tssf.false_negatives),75, ...
        'rv','filled');
    
    view(34,5);
    
    if plot_select_nuc
        scatter3(a0_for_plot(nuc_id),mRNA_for_plot(nuc_id),prot_for_plot(nuc_id),200,'ko');
        for ii = 1:nnuc_id
            text(a0_for_plot(nuc_id(ii)),mRNA_for_plot(nuc_id(ii)),prot_for_plot(nuc_id(ii)), ...
                ['   ' num2str(ii)],'fontsize',12)
        end
    end
    
    hold off;
    
    l = legend([h1,h2,h3,h4],{'low Ftz','high Ftz', ...
        'false negative','false positive'},'location','northeast');
    cur_pos = get(gcf,'position');
    set(gcf,'position',[cur_pos(1),cur_pos(2),500,400]);
    set_figure_defaults(gcf,'font')
    set(l,'box','on')
    
    xlim([0,5e5])
    zlim([0,2.5e6])
    
    
    %% -- SIMULATED VS. EMPIRICAL TRACES (ELONGATED) -- %
    tlims = [time_point_final(1),120];
    
    figure
    for ii = 1:nnuc_id
        subplot(nnuc_id,1,ii)
        plot(tlims,prot_on_thresh*[1,1], ...
            '--','color',0.7*[1 1 1],'linewidth',1.5); hold on;
        hs = plot(time_point_final,protein_trace(nuc_id(ii),:),'color',0.5*color_ftz_green);
        %         hs = [hs,plot(time_point_sim,Ptot_all(:,nuc_id(ii)),'k','linewidth',1.5)];
        hs = [hs,plot(time_point_sim,Pearly_all(:,nuc_id(ii)),'color',color_early_blue)];
        hs = [hs,plot(time_point_sim,Plate_all(:,nuc_id(ii)),'color',color_autoreg_red)];
        hs = [hs,plot(time_point_sim,Ptot_all(:,nuc_id(ii)),'color',color_ftz_green)];
        hold off;
        title(num2str(ii))
        xlim(tlims);
        ylim([0,3e6])
        ylabel('total Ftz (a.u.)')
        
        if ii == 1
            legend(hs,{'empirical','early Ftz','late Ftz','total Ftz'}, ...
                'location','northeast','box','off');   
        elseif ii == nnuc_id
            xlabel('time (min)')
        end
    end
    hold off
    cur_pos = get(gcf,'position');
    set(gcf,'position',[1070,60*nnuc_id,400,200*nnuc_id]);
    set(gcf,'units','centimeters');
    set_figure_defaults(gcf);
    
    
    %% -- ELONGATED SCHEMATIC EARLY VS. LATE PROTEIN -- %%
    sample_id = 84;
    tlims = [-20,120];
    
    figure
    plot(tlims,prot_on_thresh*[1,1], ...
        '--','color',0.7*[1 1 1],'linewidth',1.5); hold on;
    hs = plot(time_point_sim,Pearly_all(:,sample_id),'color',color_early_blue);
    hs = [hs,plot(time_point_sim,Plate_all(:,sample_id),'color',color_autoreg_red)];
    hs = [hs,plot(time_point_sim,Ptot_all(:,sample_id),'color',color_ftz_green)];
    hold off;
    xlim(tlims);
    ylim([0,3e6])
    ylabel('Ftz (a.u.)')
    xlabel('time (min)')
    
    legend(hs,'early Ftz','late Ftz','total Ftz')
    set(gcf,'position',[560 768 550 180])
    set_figure_defaults(gcf)
    
    
    %% IDENTIFY BEST PREDICTED NUCLEI
    [~,ix_sim_gastrulation] = min(abs(time_point_sim));
    
    cum_error = nansum(sqrt(((Rtot_all(1:ix_sim_gastrulation,:) - ...
        mRNA(ix_sim_start:ix_gastrulation,:)).^2 + ...
        (Ptot_all(1:ix_sim_gastrulation,:) - ...
        prot(ix_sim_start:ix_gastrulation,:)).^2)),1);
    
    err_cutoff = 7e7;
    
    figure
    hist(cum_error)
    hold on
    plot(err_cutoff*[1,1],[0,50],'r--')
    hold off
    xlabel('cumulative error')
    ylabel('count')
    set_figure_defaults(gcf);
    
    best_predicted_ix = (cum_error < err_cutoff);
    best_predicted_on = ix_on(best_predicted_ix);
    
    nbp2plt = 8;
    % find(best_predicted_ix); % choose sample indices from here
    sample_bp_ix = [2,8,15,29,48,80,110,118];
    
    Pemp_to_plot = protein_trace(sample_bp_ix,:)';
    Psim_to_plot = Ptot_all(:,sample_bp_ix);
    
    figure
    for ii_nuc = 1:nbp2plt
        subplot(2,4,ii_nuc)
        plot([time_point_final(1),0],prot_on_thresh*[1,1], ...
            '--','color',0.7*[1 1 1],'linewidth',1.5); hold on;
        plot(time_point_final,Pemp_to_plot(:,ii_nuc),'b');
        plot(time_point_sim,Psim_to_plot(:,ii_nuc),'k','linewidth',1.5);
        hold off;
        xlim([time_point_final(1),0]);
        ylim([0,3e6])
        ylabel('total Ftz (a.u.)')
        axis square
    end
    set(gcf,'position',[550 650 550 275])
    set_figure_defaults(gcf)
    
    
    sample_nbp_ix = [5,11,17,54,60,65,70,77];
    Pemp_to_plot = protein_trace(sample_nbp_ix,:)';
    Psim_to_plot = Ptot_all(:,sample_nbp_ix);
    figure
    for ii_nuc = 1:nbp2plt
        subplot(2,4,ii_nuc)
        plot([time_point_final(1),0],prot_on_thresh*[1,1], ...
            '--','color',0.7*[1 1 1],'linewidth',1.5); hold on;
        plot(time_point_final,Pemp_to_plot(:,ii_nuc),'b');
        plot(time_point_sim,Psim_to_plot(:,ii_nuc),'k','linewidth',1.5);
        hold off;
        xlim([time_point_final(1),0]);
        ylim([0,3e6])
        ylabel('total Ftz (a.u.)')
        axis square
    end
    set(gcf,'position',[550 650 550 275])
    set_figure_defaults(gcf)
    
    
    %% FRACTION CORRECTLY CLASSIFIED NUCLEI THAT MAINTAIN FATE
    use_for_fate_analysis = best_predicted_ix & tssf.correctly_identified;
    
    maintain_fate = use_for_fate_analysis & ~ix_disagree;
    maintain_on_fate = ix_on & maintain_fate;
    maintain_off_fate = ~ix_on & maintain_fate;
    disp(['Correctly predicted that maintain fate: ' ...
        num2str(100*sum(maintain_fate)/sum(use_for_fate_analysis)) ...
        '% (' num2str(sum(maintain_fate)) ' of ' ...
        num2str(sum(use_for_fate_analysis)) ...
        ' nuclei)']);
    disp([' - ' ...
        num2str(100*sum(maintain_on_fate)/sum(ix_on & use_for_fate_analysis)) ...
        '% of on nuclei (' num2str(sum(maintain_on_fate)) ' of ' ...
        num2str(sum(ix_on & use_for_fate_analysis)) ')']);
    disp([' - ' ...
        num2str(100*sum(maintain_off_fate)/sum(~ix_on & use_for_fate_analysis)) ...
        '% of off nuclei (' num2str(sum(maintain_off_fate)) ' of ' ...
        num2str(sum(~ix_on & use_for_fate_analysis)) ')']);
    
    
    %% VISUALIZE SIMULATED VS. EMPIRICAL EMBRYOS
    
    nucsize = 50;
    
    xcoord_valid = xcoord(valid_ix);
    ycoord_valid = ycoord(valid_ix);
    
    false_negative_no_autoreg = pred_off_thresh_no_autoreg(ix_sim_gastrulation,:) & ix_on;
    false_positive_no_autoreg = pred_on_thresh_no_autoreg(ix_sim_gastrulation,:) & ~ix_on;
    
    times_to_plot = [-20,-10,0];
    ntpts = numel(times_to_plot);
    
    emb_to_plot = unique(embryoid);
    
    show_false_ids = true;
    
    % normalize all to s.s. protein value so comparable
    [~,ssP] = fi.get_steady_states(inf);
    norm_factor = ssP(end);
    emp_prot_normalized = prot / norm_factor;
    sim_prot_normalized = Ptot_all / norm_factor;
    sim_prot_no_autoreg_normalized = Pearly_all / norm_factor;
    for ii_emb = 1:numel(emb_to_plot)
        cur_embryoid = emb_to_plot(ii_emb);
        if ~ismember(cur_embryoid,embryoid_to_ignore)
            figure('name',['embryo ' num2str(cur_embryoid)]);
            cur_emb_ix = (embryoid == cur_embryoid)';
            cur_ix = (valid_ix & cur_emb_ix);
            for ii_t = 1:ntpts
                [~,empix2plt] = min(abs(time_point_final - times_to_plot(ii_t)));
                [~,simix2plt] = min(abs(time_point_sim - times_to_plot(ii_t)));
                
                subplot(3,ntpts,ii_t)
                scatter(rot_emb(emb_to_plot(ii_emb))*xcoord(cur_ix),rot_emb(emb_to_plot(ii_emb))*ycoord(cur_ix),nucsize, ...
                    emp_prot_normalized(empix2plt, ...
                    cur_ix(valid_ix))'*[0 1 0],'filled');
                if show_false_ids && ii_t == ntpts
                    cur_ix_do_not_maintain_fate = cur_ix(valid_ix) & use_for_fate_analysis & ~maintain_fate;
                    hold on
                    scatter(rot_emb(emb_to_plot(ii_emb))*xcoord_valid(cur_ix_do_not_maintain_fate), ...
                        rot_emb(emb_to_plot(ii_emb))*ycoord_valid(cur_ix_do_not_maintain_fate),50,'rx')
                    hold off
                end
                title(['t = ' num2str(times_to_plot(ii_t))])
                axis equal
                axis off
                
                subplot(3,ntpts,ntpts+ii_t)
                scatter(rot_emb(emb_to_plot(ii_emb))*xcoord(cur_ix),rot_emb(emb_to_plot(ii_emb))*ycoord(cur_ix),nucsize, ...
                    sim_prot_normalized(simix2plt, ...
                    cur_ix(valid_ix))'*[0 1 0],'filled');
                if show_false_ids && ii_t == ntpts
                    hold on
                    scatter(rot_emb(emb_to_plot(ii_emb))*xcoord_valid(cur_ix(valid_ix) & tssf.false_negatives), ...
                        rot_emb(emb_to_plot(ii_emb))*ycoord_valid(cur_ix(valid_ix) & tssf.false_negatives),50,'r^')
                    scatter(rot_emb(emb_to_plot(ii_emb))*xcoord_valid(cur_ix(valid_ix) & tssf.false_positives), ...
                        rot_emb(emb_to_plot(ii_emb))*ycoord_valid(cur_ix(valid_ix) & tssf.false_positives),50,'rv')
                    hold off
                end
                axis equal
                axis off
                
                subplot(3,ntpts,2*ntpts+ii_t)
                scatter(rot_emb(emb_to_plot(ii_emb))*xcoord(cur_ix),rot_emb(emb_to_plot(ii_emb))*ycoord(cur_ix),nucsize, ...
                    sim_prot_no_autoreg_normalized(simix2plt, ...
                    cur_ix(valid_ix))'*[0 1 0],'filled');
                if show_false_ids && ii_t == ntpts
                    hold on
                    scatter(rot_emb(emb_to_plot(ii_emb))*xcoord_valid(cur_ix(valid_ix) & false_negative_no_autoreg), ...
                        rot_emb(emb_to_plot(ii_emb))*ycoord_valid(cur_ix(valid_ix) & false_negative_no_autoreg),50,'r^')
                    scatter(rot_emb(emb_to_plot(ii_emb))*xcoord_valid(cur_ix(valid_ix) & false_positive_no_autoreg), ...
                        rot_emb(emb_to_plot(ii_emb))*ycoord_valid(cur_ix(valid_ix) & false_positive_no_autoreg),50,'rv')
                    hold off
                end
                axis equal
                axis off
            end
            set(gcf,'Position',[60 0 580 430]);
        end
    end
    
    
    %% ESTIMATED CONVERGENCE RATE TO S.S.
    fprime = @(x) ((fi.a*fi.n*x.^(fi.n-1))*fi.K^fi.n)./((fi.K^fi.n + x.^fi.n).^2);
    lambda_pos = (-(fi.gammaP + fi.gammaR) + sqrt((fi.gammaP + fi.gammaR)^2 - ...
        4*(fi.gammaR*fi.gammaP - fi.alpha*fi.c*fprime(ssP))))/2;
    lambda_neg = (-(fi.gammaP + fi.gammaR) - sqrt((fi.gammaP + fi.gammaR)^2 - ...
        4*(fi.gammaR*fi.gammaP - fi.alpha*fi.c*fprime(ssP))))/2;
    
    con_rate = min([abs(lambda_pos(:)),fi.b*[1;1;1],fi.gammaR*[1;1;1]],2);
    con_halflife = log(2)./con_rate;   % min
    disp(['Convergence rate "half-life" for low state = ' ...
        num2str(con_halflife(1)) ' min'])
    disp(['Convergence rate "half-life" for high state = ' ...
        num2str(con_halflife(3)) ' min'])
    
    
    %% ESTIMATED ACCURACY WHEN AUTOREGULATORY ELEMENT DELAYED TURNING ON
    % only consider correctly identified nuclei that are fated on (predictions
    % for off will never change through early mRNA decay)
    fated_on = use_for_fate_analysis & ssf.predicted_on;
    
    % determine when crosses separatrix
    timepts_above_surface = false(size(Pearly_all));
    time2xing = nan(numel(fated_on),1);
    for ii = 1:numel(fated_on)
        if fated_on(ii)
            timepts_above_surface(:,ii) = ssf.above_surface( ...
                aearly_all(:,ii),Rearly_all(:,ii),Pearly_all(:,ii));
            xix = find(~timepts_above_surface(:,ii),1);
            if ~isempty(xix)
                time2xing(ii) = tarr(xix);
            end
        end
    end
    
    
    cols_by_time = (1-(time2xing - nanmin(time2xing))/range(time2xing))*[1 0 0];
    cols_by_time(isnan(time2xing)' & fated_on,:) = 0.7;
    
    figure
    ss = surf(ssf.Xq,ssf.Yq,ssf.Zq);
    shading interp;
    set(ss,'facealpha',0.5,'facecolor','k');
    xlabel('a_0'); ylabel('r_0'); zlabel('p_0');
    ylim([0,3e6]);
    
    hold on;
    plot3(aearly_all(:,fated_on), ...
        Rearly_all(:,fated_on), ...
        Pearly_all(:,fated_on),'--','linewidth',0.5, ...
        'color',0.7*[1 1 1]);
    
    for ii = 1:numel(fated_on)
        if fated_on(ii)
            plot3(aearly_all(~timepts_above_surface(:,ii),ii), ...
                Rearly_all(~timepts_above_surface(:,ii),ii), ...
                Pearly_all(~timepts_above_surface(:,ii),ii), ...
                'color',cols_by_time(ii,:));
        end
    end
    
    scatter3(a0(fated_on), ...
        mRNA_for_plot(fated_on), ...
        prot_for_plot(fated_on),75, ...
        cols_by_time(fated_on,:),'filled');
    
    hold off;
    
    view(-5,40);
    xlim([0,5e5])
    
    correctly_predicted_no_autoreg = sum(isnan(time2xing)' & fated_on);
    
    disp(['Autoregulatory element can be delayed by ' ...
        num2str(nanmin(time2xing)) ' min without changing any predictions'])
    disp(['   ' num2str(sum(fated_on)) ' nuclei total, ' ...
        num2str(correctly_predicted_no_autoreg) ' (' ...
        num2str(100*correctly_predicted_no_autoreg/sum(fated_on)) ...
        '%) do not change fate for delay up to ' num2str(max(tarr)) ' min']);
    
    
    [est_cdf,est_cdf_eval_pts] = ecdf(time2xing);
    proportion_change_fate = sum(~isnan(time2xing))/sum(fated_on);
    set_figure_defaults(gcf,'font');
    
    xlabel('r_0')
    ylabel('R_0')
    zlabel('P_0')
    
    figure
    plot(est_cdf_eval_pts + time_autoreg_start,1-est_cdf*proportion_change_fate, ...
        'color',color_latestart_gray)
    xlabel({'t_{on} (min)'})
    ylabel({'fraction', 'fated on'})
    axis square
    xlim([0,85]+time_autoreg_start);
    ylim([0,1])
    box off
    set_figure_defaults(gcf)
    set(gcf,'position',[550 800 250 160])
    
    
    %% ESTIMATED TIME EARLY ELEMENT CAN CEASE PRODUCTION
    % find time when total protein falls within basin of attraction to high
    % steady state for autoregulatory element
    
    % identify separatrix
    R_set = [0,3.27e6];
    P_set = [2.42e6,0];
    fi.set('a0',0,'r0',0,'p0',0);
    Rsep_ne = [];
    Psep_ne = [];
    Rsep_nw = [];
    Psep_nw = [];
    tsep = 0:0.2:200;
    for ii = 1:length(R_set)
        [M,P,t] = fi.simulate([R_set(ii),P_set(ii)],tsep, ...
            'no_plot','no_save','off');
        switch ii
            case 1
                Rsep_nw = [Rsep_nw;M];
                Psep_nw = [Psep_nw;P];
            case 2
                Rsep_ne = [Rsep_ne;M];
                Psep_ne = [Psep_ne;P];
        end
    end
    
    % find where first approach unstable s.s. and remove all points thereafter
    [ssR,ssP] = fi.get_steady_states(inf);
    
    err_ne = sqrt((ssR(2) - Rsep_ne).^2 + (ssP(2) - Psep_ne).^2);
    [~,ix_min] = min(err_ne);
    Rsep_ne(ix_min:end) = [];
    Psep_ne(ix_min:end) = [];
    
    err_nw = sqrt((ssR(2) - Rsep_nw).^2 + (ssP(2) - Psep_nw).^2);
    [~,ix_min] = min(err_nw);
    Rsep_nw(ix_min:end) = [];
    Psep_nw(ix_min:end) = [];
    
    Rsep = [Rsep_ne;flipud(Rsep_nw)];
    Psep = [Psep_ne;flipud(Psep_nw)];
    
    % smooth boundary
    smspan = 10;
    Rsep = smooth(Rsep,smspan);
    Psep = smooth(Psep,smspan);
    
    Rsep = [flipud(Rsep);max(max(Rtot_all))];
    Psep = [flipud(Psep);0];
    
    % determine whether total mRNA and protein fall within basin of attraction
    % for high Ftz
    in_boa = (Ptot_all > interp1(Rsep,Psep,Rtot_all));
    
    figure
    plot(time_point_sim,sum(in_boa(:,maintain_on_fate)/sum(maintain_on_fate),2), ...
        'color',color_shutdown_gold);
    set_figure_defaults(gcf)
    xlim([time_point_sim(1),65]);
    ylim([0,1])
    axis square
    xlabel('t_{off} (min)')
    ylabel({'fraction', 'fated on'})
    set(gcf,'position',[550 800 250 160])
 
    
    %% ESTIMATED GAP BETWEEN EARLY SHUTOFF AND AUTOREGULATORY START
    if estimate_gap
        %% GENERATE HEATMAP
        nfatedon = sum(fated_on);
        nuc2sim = find(fated_on);
        
        ixs_delay = 1:10:ix_mirrorautoreg;
        times_early_shutdown = time_point_sim(ixs_delay);
        times_late_start = time_point_sim(ixs_delay);
        times_delay = times_late_start - times_early_shutdown';
        ee_before_ar = times_delay >= 0;
        ee_after_ar = ~ee_before_ar;
        
        fi_delay = ftz_instance();
        
        % rows: late start, cols: shutdown
        R_at_t = nan(numel(times_late_start),numel(times_early_shutdown),nfatedon);
        P_at_t = nan(numel(times_late_start),numel(times_early_shutdown),nfatedon);
        
        for ii_nuc = 1:size(R_at_t,3)
            for ii_sd = 1:numel(times_early_shutdown)
                cur_delays = times_late_start - times_early_shutdown(ii_sd);
                ix_cur_ee_first = (cur_delays >= 0);
                
                fi_delay.set('a0',0,'r0',Rearly_all(ixs_delay(ii_sd),nuc2sim(ii_nuc)), ...
                    'p0',Pearly_all(ixs_delay(ii_sd),nuc2sim(ii_nuc)));
                R_at_t(ix_cur_ee_first,ii_sd,ii_nuc) = fi_delay.r(cur_delays(ix_cur_ee_first));
                P_at_t(ix_cur_ee_first,ii_sd,ii_nuc) = fi_delay.p(cur_delays(ix_cur_ee_first));
            end
        end
        
        for ii_nuc = 1:size(R_at_t,3)
            for ii_ls = 1:numel(times_late_start)
                for ii_sd = 1:numel(times_early_shutdown)
                    cur_delay = times_late_start(ii_ls) - times_early_shutdown(ii_sd);
                    if cur_delay < 0
                        % early element turns off after autoregulation turns on
                        
                        % forward simulate from autoreg on to shutdown
                        fi_delay.set('a0',aearly_all(ixs_delay(ii_ls),nuc2sim(ii_nuc)), ...
                            'r0',Rearly_all(ixs_delay(ii_ls),nuc2sim(ii_nuc)), ...
                            'p0',Pearly_all(ixs_delay(ii_ls),nuc2sim(ii_nuc)));
                        [R,P] = fi_delay.simulate([0,0],[0,-cur_delay],'no_plot','no_save');
                        R_at_t(ii_ls,ii_sd,ii_nuc) = R(end) + fi_delay.r(-cur_delay);
                        P_at_t(ii_ls,ii_sd,ii_nuc) = P(end) + fi_delay.p(-cur_delay);
                    end
                end
            end
        end
        
        delay_in_boa = (P_at_t > interp1(Rsep,Psep,R_at_t));
        delay_frac = squeeze(sum(delay_in_boa,3)/nfatedon);
        
        figure
        % x to column indices, y to row indices
        surf(times_early_shutdown,times_late_start,delay_frac, ...
            'EdgeColor','none');
        view([0,90]);
        colorbar('northoutside'); colormap jet;
        xticks(-20:5:20)
        yticks(-20:5:20)
        xlim([min(times_early_shutdown),max(times_early_shutdown)]);
        ylim([min(times_late_start),max(times_late_start)]);
        xlabel('time early element production ceases (min)')
        ylabel('time autoregulatory element becomes responsive (min)')
        axis square
        set_figure_defaults(gcf);
        
        
        %% PLOT HEATMAP
        figure
        % x to column indices, y to row indices
        surf(times_early_shutdown,times_late_start,delay_frac, ...
            'EdgeColor','none');
        
        view([0,90]);
        caxis([0,1])
        colorbar('northoutside','Ticks',[0 0.5 1]); colormap gray;
        xticks(-20:5:20)
        yticks(-20:5:20)
        xlim([min(times_early_shutdown),max(times_early_shutdown)]);
        ylim([min(times_late_start),max(times_late_start)]);
        xlabel('time early element production ceases (min)')
        ylabel({'time autoregulatory element becomes', 'responsive (min)'})
        axis square
        set_figure_defaults(gcf);
        
        title('fraction nuclei adopting correct fate')
    end
end


%% COMMITMENT WINDOW
% only consider cases where autoregulation start precedes early element
% shutdown (since most nuclei with the opposite adopt the low fate)
if estimate_window
    nfatedon = sum(fated_on);
    nuc2sim = find(fated_on);
    
    [~,ix_delay_end] = min(abs(time_point_sim - 40));
    ixs_sd_delay = 1:5:ix_delay_end;
    ixs_ls_delay = 1:5:ix_sim_gastrulation;
    times_early_shutdown = time_point_sim(ixs_sd_delay);
    times_late_start = time_point_sim(ixs_ls_delay);
    times_delay = times_late_start - times_early_shutdown';
    ee_before_ar = times_delay >= 0;
    ee_after_ar = ~ee_before_ar;
    
    fi_delay = ftz_instance();
    
    % rows: late start, cols: shutdown
    R_at_t = nan(numel(times_late_start),numel(times_early_shutdown),nfatedon);
    P_at_t = nan(numel(times_late_start),numel(times_early_shutdown),nfatedon);
    
    for ii_nuc = 1:size(R_at_t,3)
        for ii_ls = 1:numel(times_late_start)
            for ii_sd = 1:numel(times_early_shutdown)
                cur_delay = times_late_start(ii_ls) - times_early_shutdown(ii_sd);
                if cur_delay < 0
                    % early element turns off after autoregulation turns on
                    
                    % forward simulate from autoreg on to shutdown
                    fi_delay.set('a0',aearly_all(ixs_ls_delay(ii_ls),nuc2sim(ii_nuc)), ...
                        'r0',Rearly_all(ixs_ls_delay(ii_ls),nuc2sim(ii_nuc)), ...
                        'p0',Pearly_all(ixs_ls_delay(ii_ls),nuc2sim(ii_nuc)));
                    [R,P] = fi_delay.simulate([0,0],[0,-cur_delay],'no_plot','no_save');
                    R_at_t(ii_ls,ii_sd,ii_nuc) = R(end) + fi_delay.r(-cur_delay);
                    P_at_t(ii_ls,ii_sd,ii_nuc) = P(end) + fi_delay.p(-cur_delay);
                end
            end
        end
    end
    
    delay_in_boa = (P_at_t > interp1(Rsep,Psep,R_at_t));
    delay_frac = squeeze(sum(delay_in_boa,3)/nfatedon);
    
    
    %% PLOT EXAMPLE TRACES
    sd_set = [10,20,30,40,50];
    ls_set = [5,5,15,5,21];
    
    cur_nuc_mf = 21;
    temp_ixs = find(fated_on);
    cur_nuc = temp_ixs(cur_nuc_mf);
    
    ymax = 3e6;
    
    figure
    xlim([0,50] + time_point_sim(1))
    ylim([0,ymax])
    axis square
    xlabel('time (min)')
    ylabel('Ftz')
    
    for ii_set = 1:numel(sd_set)
        ii_cur_sd = sd_set(ii_set);
        
        ii_cur_ls = ls_set(ii_set);
        
        time_delay_sim = time_point_sim(1:end);
        [~,cur_sd] = min(abs(time_delay_sim - times_early_shutdown(ii_cur_sd)));
        [~,cur_ls] = min(abs(time_delay_sim - times_late_start(ii_cur_ls)));
        time_delay_sim = time_delay_sim - time_delay_sim(1);
        
        Rdelay = nan(numel(time_delay_sim),1);
        Pdelay = nan(numel(time_delay_sim),1);
        cur_delay = times_late_start(ii_cur_ls) - times_early_shutdown(ii_cur_sd);
        if cur_delay < 0
            fi_delay.set('a0',a(ix_sim_start,cur_nuc), ...
                'r0',mRNAearly(ix_sim_start,cur_nuc), ...
                'p0',Pearly(ix_sim_start,cur_nuc));
            
            Rdelay(1:cur_ls) = fi_delay.r(time_delay_sim(1:cur_ls));
            Pdelay(1:cur_ls) = fi_delay.p(time_delay_sim(1:cur_ls));
            
            fi_delay.set('a0',fi_delay.c1*exp(-fi_delay.b*time_delay_sim(cur_ls)), ...
                'r0',Rdelay(cur_ls),'p0',Pdelay(cur_ls));
            [R,P] = fi_delay.simulate([0,0],time_delay_sim(cur_ls:cur_sd) - time_delay_sim(cur_ls),'no_plot','no_save');
            
            Rdelay(cur_ls:cur_sd) = R + fi_delay.r(time_delay_sim(cur_ls:cur_sd) - time_delay_sim(cur_ls))';
            Pdelay(cur_ls:cur_sd) = P + fi_delay.p(time_delay_sim(cur_ls:cur_sd) - time_delay_sim(cur_ls))';
            
            fi_delay.set('a0',0,'r0',0,'p0',0);
            
            if cur_sd ~= length(time_delay_sim)
                [R,P] = fi_delay.simulate([Rdelay(cur_sd),Pdelay(cur_sd)], ...
                    time_delay_sim(cur_sd:end) - time_delay_sim(cur_sd),'no_plot','no_save');
                Rdelay(cur_sd:end) = R + fi_delay.r(time_delay_sim(cur_sd:end) - time_delay_sim(cur_sd))';
                Pdelay(cur_sd:end) = P + fi_delay.p(time_delay_sim(cur_sd:end) - time_delay_sim(cur_sd))';
            end
        else
            % early element turns off before autoregulation turns on
            fi_delay.set('a0',a(ix_sim_start,cur_nuc), ...
                'r0',mRNAearly(ix_sim_start,cur_nuc), ...
                'p0',Pearly(ix_sim_start,cur_nuc));
            
            Rdelay(1:cur_sd) = fi_delay.r(time_delay_sim(1:cur_sd));
            Pdelay(1:cur_sd) = fi_delay.p(time_delay_sim(1:cur_sd));
            
            fi_delay.set('a0',0,'r0',Rdelay(cur_sd),'p0',Pdelay(cur_sd));
            
            Rdelay(cur_sd:cur_ls) = fi_delay.r(time_delay_sim(cur_sd:cur_ls) - time_delay_sim(cur_sd));
            Pdelay(cur_sd:cur_ls) = fi_delay.p(time_delay_sim(cur_sd:cur_ls) - time_delay_sim(cur_sd));
            
            fi_delay.set('r0',Rdelay(cur_ls),'p0',Pdelay(cur_ls));
            
            [R,P] = fi_delay.simulate([0,0],time_delay_sim(cur_ls:end) - time_delay_sim(cur_ls),'no_plot','no_save');
            Rdelay(cur_ls:end) = R + fi_delay.r(time_delay_sim(cur_ls:end) - time_delay_sim(cur_ls))';
            Pdelay(cur_ls:end) = P + fi_delay.p(time_delay_sim(cur_ls:end) - time_delay_sim(cur_ls))';
        end
        
        subplot(1,numel(sd_set),ii_set)
        plot(time_delay_sim + time_point_sim(1),Pdelay,'color',color_ftz_green);
        hold on;
        plot(times_early_shutdown(ii_cur_sd)*[1,1],[0,ymax],'color',color_shutdown_gold);
        plot(times_late_start(ii_cur_ls)*[1,1],[0,ymax],'color',color_latestart_gray);
        hold off;
        xlim([time_delay_sim(1),time_delay_sim(end)] + time_point_sim(1))
        ylim([0,ymax])
        xlabel('time (min)')
        ylabel('Ftz')
    end
    set_figure_defaults(gcf);
    set(gcf,'position',[550 800 1080 145]);
    
    
    %% CONTOUR PLOT
    % error tolerance to define window
    errtols = [0.05,0.15,0.25,0.5];
    
    comm_win = nan(numel(times_late_start),numel(errtols));
    for ii_tol = 1:numel(errtols)
        tol = errtols(ii_tol);
        
        % autoreg start time vs. toff - ton to get no more than tol errors
        for ii_ls = 1:numel(times_late_start)
            valix = (1 - delay_frac(ii_ls,:)) <= tol;
            if any(valix)
                comm_win(ii_ls,ii_tol) = min(times_early_shutdown(valix)) - ...
                    times_late_start(ii_ls);
            end
        end
    end
    
    figure
    plot(times_late_start,comm_win)
    xlim([times_late_start(1),0])
    ylim([0,35])
    xlabel('autoregulation start')
    ylabel('commitment window')
    
    set(gcf,'position',[2300 750 150 150])
end


%% ESTIMATES FOR RNA AND PROTEIN NOISE

if get_noise_params
    dt = time_point_final(2)-time_point_final(1);
    
    dRlate = dRNA_trace(valid_ix,:);
    mRNA_late_est = cumsum(dRlate*dt,2);
    mRNA_late_est(mRNA_late_est < 0) = NaN;
    mRNA_late_est(:,time_point_final < time_autoreg_start) = NaN;
    
    dRlate_dt_term = real(fi.c*fi.f(0, ...
        max(protein_trace,fi.p(inf)))) - fi.gammaR*mRNA_late_est;
    dRlate_noise_term = dRlate - dRlate_dt_term*dt;
    
    mRNA_nbin = 8;
    quant_edges = unique(quantile(mRNA_late_est(:),mRNA_nbin));
    mRNA_bin = discretize(mRNA_late_est,quant_edges);
    
    muR = nan(mRNA_nbin,1);
    pdR = prob.NormalDistribution.empty(mRNA_nbin,0);
    invalid_ix = false(mRNA_nbin,1);
    figure
    for jj_bin = 1:mRNA_nbin
        muR(jj_bin) = nanmean(mRNA_late_est(mRNA_bin == jj_bin));
        cur_vals = dRlate_noise_term(mRNA_bin == jj_bin)/sqrt(dt);
        histogram(cur_vals + muR(jj_bin), ...
            'DisplayStyle','stairs','normalization','pdf');
        hold on;
        
        nanix = isnan(cur_vals);
        cur_vals(nanix) = [];
        if ~isempty(cur_vals)
            pdR(jj_bin) = fitdist(cur_vals + muR(jj_bin),'normal');
        else
            invalid_ix(jj_bin) = true;
        end
    end
    hold off;
    
    muR(invalid_ix) = [];
    
    xlabel('late mRNA (A.U.)')
    ylabel('pdf')
    title({['Distributions of stochastic contribution to noise binned'], ...
        ['and centered at mean late mRNA value per bin']});
    set_figure_defaults(gcf);
    
    
    mu_to_fit = [pdR.mu]';
    var_to_fit = [pdR.sigma]'.^2;
    
    muRvar_fit = fit(mu_to_fit,var_to_fit,fittype('a*x + b'), ...
        'Lower',[0 0],'Upper',[inf inf],'StartPoint',[1e5 0]);
    
    figure
    scatter(mu_to_fit,var_to_fit); hold on;
    spt = linspace(0,1e6,100);
    plot(spt,muRvar_fit(spt)); hold off;
    l = legend('data','fit');
    set(l,'location','southeast')
    xlabel('mean late mRNA level')
    ylabel('variance')
    title({'Variance of process noise as function of ' 'mRNA concentration'})
    
    set_figure_defaults(gcf)
    
    early_mRNA_determ = nan(nvalid_nuc,numel(time_point_final));
    tt = time_point_final;
    tt(time_point_final < time_autoreg_start) = NaN;
    tt = tt - time_autoreg_start;
    for ii = 1:nvalid_nuc
        fi.set('a0',a(ix_sim_start,ii),'r0',mRNAearly(ix_sim_start,ii),'p0',Pearly(ix_sim_start,ii));
        cur_early_mRNA = (1/(fi.gammaR-fi.b))*fi.c1*exp(-fi.b*tt) + ...
            fi.c2*exp(-fi.gammaR*tt);
        early_mRNA_determ(ii,:) = cur_early_mRNA;
    end
    
    
    %% PROTEIN NOISE
    valid_prot = protein_trace(:,1:end-1);
    valid_prot(valid_prot < 0) = NaN;
    dPtot = diff(protein_trace,1,2);
    dP_dt_term = fi.alpha*(early_mRNA_determ(:,1:end-1)+mRNA_late_est(:,1:end-1)) - ...
        fi.gammaP*valid_prot;
    dP_noise_term = (dPtot - dP_dt_term*dt);
    
    prot_nbin = 6;
    quant_edges = unique(quantile(valid_prot(:),prot_nbin));
    prot_bin = discretize(valid_prot,quant_edges);
    
    mu = nan(prot_nbin,1);
    pd = prob.NormalDistribution.empty(0,prot_nbin);
    figure
    for jj_bin = 1:prot_nbin
        mu(jj_bin) = mean(valid_prot(prot_bin == jj_bin));
        cur_vals = dP_noise_term(prot_bin == jj_bin)/sqrt(dt);
        histogram(cur_vals + mu(jj_bin), ...
            'DisplayStyle','stairs','normalization','pdf');
        hold on;
        
        nanix = isnan(cur_vals);
        cur_vals(nanix) = [];
        if ~isempty(cur_vals)
            pd(jj_bin) = fitdist(cur_vals + mu(jj_bin),'normal');
        end
        
    end
    hold off;
    
    xlabel('protein (A.U.)')
    ylabel('pdf')
    title({['Distributions of stochastic contribution to noise binned'], ...
        ['and centered at mean protein value per bin']});
    set_figure_defaults(gcf);
    
    
    muP_to_fit = [pd.mu]';
    varP_to_fit = [pd.sigma]'.^2;
    
    muPvar_fit = fit(muP_to_fit,varP_to_fit,fittype('a*x+b'),'Lower',[0,0], ...
        'StartPoint',[0 1e5]);
    
    figure
    scatter(muP_to_fit,varP_to_fit); hold on;
    spt = linspace(0,2e6,100);
    plot(spt,muPvar_fit(spt)); hold off;
    l = legend('data','fit');
    set(l,'location','southeast')
    xlabel('mean protein level')
    ylabel('variance')
    title({'Variance of process noise as function of ' 'protein concentration'})
    
    set_figure_defaults(gcf)
    
    if save_noise_params
        save('muRvar_fit','muRvar_fit');
        save('muPvar_fit','muPvar_fit');
    end
    
    
    sample_pts = linspace(0,2e6,200);
    
    figure
    for jj_bin = 1:numel(pd)
        if ~isempty(pd(jj_bin))
            sample_pdf = (1/(pd(jj_bin).sigma*sqrt(2*pi)))* ...
                exp(-0.5*((sample_pts-pd(jj_bin).mu)/pd(jj_bin).sigma).^2);
            plot(sample_pts,sample_pdf);
            hold on;
        end
    end
    hold off;
    xlabel('protein (au)')
    ylabel('pdf')
    set_figure_defaults(gcf)
end


%% RUN STOCHASTIC SIMULATIONS
if stochastic_sims
    tssf_stoc = switching_separatrix('trans_switching_separatrix_a0-r0-p0.mat');
    tssf_stoc.prediction_accuracy(a(ix_autoreg_start,:),mRNAearly(ix_autoreg_start,:), ...
        Pearly(ix_autoreg_start,:),ix_on,valid_ix,embryoid);
    
    %% STOCHASTIC SIMULATION BY SAMPLING FROM INSIDE PT CLOUD
    
    % only valid for boundary nuclei
    load('muRvar_fit.mat');
    load('muPvar_fit.mat');
    muRsig_fit = @(x) sqrt(muRvar_fit(x));
    muPsig_fit = @(x) sqrt(muPvar_fit(x));
    
    % make sample "experiments"
    nnuc = round(numel(a(ix_autoreg_start,:))/2);
    emp_above_surf = tssf_stoc.above_surface(a(ix_autoreg_start,:), ...
        mRNAearly(ix_autoreg_start,:),Pearly(ix_autoreg_start,:));
    nexps = 100;
    
    uvw_temp = sample_inside_pt_cloud([a(ix_autoreg_start,:)', ...
        mRNAearly(ix_autoreg_start,:)',Pearly(ix_autoreg_start,:)'],2*nexps*nnuc);
    to_draw_from_above = tssf_stoc.above_surface(uvw_temp(:,1),uvw_temp(:,2),uvw_temp(:,3));
    to_draw_from_above = uvw_temp(to_draw_from_above,:);
    uvw_above = uvw_temp(randperm(size(uvw_temp,1),nexps*nnuc),:);
    
    uvw_below = sample_inside_pt_cloud([a(ix_autoreg_start,:)', ...
        mRNAearly(ix_autoreg_start,:)',Pearly(ix_autoreg_start,:)'],nexps*nnuc);
    
    uvw = [uvw_above;uvw_below];
    
    
    tmax = -time_autoreg_start;
    tarr = 0:dt:tmax;
    ntpts = numel(tarr);
    
    nsims = size(uvw,1);
    
    U = nan(ntpts,2,nsims);
    Ptot = nan(ntpts,nsims);
    for jj_sim = 1:nsims
        disp(['sample run ' num2str(jj_sim) ' of ' num2str(nsims)]);
        
        fi.set('a0',uvw(jj_sim,1),'r0',uvw(jj_sim,2),'p0',uvw(jj_sim,3));
        U(:,:,jj_sim) = fi.sde([0;uvw(jj_sim,3)],tarr,muRsig_fit,muPsig_fit);
        Ptot(:,jj_sim) = U(:,2,jj_sim);
    end
    
    disp(['Simulations ending in NaN = ' num2str(sum(isnan(Ptot(end,:)))) ...
        ' out of ' num2str(nsims)]);
    
    
    %% SAMPLE STOCHASTIC TRACES
    on_ixs = find(Ptot(end,:) > fi.prot_on_thresh);
    off_ixs = find(Ptot(end,:) < fi.prot_on_thresh);
    sample_ids = [on_ixs([1,2]),off_ixs([3,5])];
    
    figure
    subplot(2,1,1)
    plot(tarr+time_autoreg_start,squeeze(U(:,1,sample_ids)));
    xlabel('time before gastrulation (min)')
    ylabel('late mRNA (au)')
    
    subplot(2,1,2)
    plot(tarr+time_autoreg_start,Ptot(:,sample_ids));
    hold on;
    plot([min(tarr),max(tarr)]+time_autoreg_start,fi.prot_on_thresh*[1,1],'--','color',0.7*[1 1 1]);
    hold off;
    ylim([0,2e6])
    xlabel('time before gastrulation (min)')
    ylabel('total protein (au)')
    
    set_figure_defaults(gcf)

    
    % group into "experiments"
    false_pos_rate = nan(nexps,1);
    false_neg_rate = nan(nexps,1);
    tot_error_rate = nan(nexps,1);
    for ii_exp = 1:nexps
        start_ix = nnuc*(ii_exp-1);
        cur_ix = [(nnuc*(ii_exp-1)+1):(nnuc*(ii_exp-1)+nnuc), ...
            nexps*nnuc+((nnuc*(ii_exp-1)+1):(nnuc*(ii_exp-1)+nnuc))];
        
        cur_above_surface = tssf_stoc.above_surface(uvw(cur_ix,1), ...
            uvw(cur_ix,2),uvw(cur_ix,3));
        
        false_pos_rate(ii_exp) = sum(cur_above_surface & (Ptot(end,cur_ix) < prot_on_thresh)) / (2*nnuc);
        false_neg_rate(ii_exp) = sum(~cur_above_surface & (Ptot(end,cur_ix) >= prot_on_thresh)) / (2*nnuc);
        tot_error_rate(ii_exp) = false_pos_rate(ii_exp) + false_neg_rate(ii_exp);
    end
    
    
    [sim_exp_err_rate_ecdf,xx] = ecdf(tot_error_rate);
    [sim_exp_fneg_rate_ecdf,xxnr] = ecdf(false_neg_rate);
    [sim_exp_fpos_rate_ecdf,xxpr] = ecdf(false_pos_rate);

    nuc_of_interest = (embryoid(valid_ix) ~= 3);
    emp_false_neg_rate = sum(tssf_stoc.false_negatives(nuc_of_interest)) / sum(nuc_of_interest);
    emp_false_pos_rate = sum(tssf_stoc.false_positives(nuc_of_interest)) / sum(nuc_of_interest);
    emp_error_rate = sum(tssf_stoc.false_negatives(nuc_of_interest) | tssf_stoc.false_positives(nuc_of_interest)) / sum(nuc_of_interest);

    % the first value returned is for convenience in stair plotting so we
    % ignore
    plot(xx(2:end),sim_exp_err_rate_ecdf(2:end),'k'); hold on;
    plot(xxnr(2:end),sim_exp_fneg_rate_ecdf(2:end),'r');
    plot(xxpr(2:end),sim_exp_fpos_rate_ecdf(2:end),'g');
    
    ax = gca;
    ymax = ax.YLim(2);
    
    plot(emp_error_rate*[1,1],[0 ymax],'k--');
    plot(emp_false_neg_rate*[1,1],[0 ymax],'r--');
    plot(emp_false_pos_rate*[1,1],[0 ymax],'g--');
    hold off;
    
    xlabel('rate')
    ylabel('cumulative probability')
    
    legend({'total errors','false negatives','false positives'}, ...
        'location','southeast')
    
    set_figure_defaults(gcf)
    
    
    %% STOCHASTIC SWITCHING
    nsims_sw = 100;
    tmax = 420;	% 7 hours from gastrulation into germ band retraction (stage 12-13)
    tarr = 0:0.2:tmax;
    ntpts = numel(tarr);
    
    [ssR,ssP] = fi.get_steady_states(inf);
    
    fi.set('a0',0,'r0',0,'p0',0);
    U_sw_lh = nan(ntpts,2,nsims_sw);
    Ptot_sw_lh = nan(ntpts,nsims_sw);
    U_sw_hl = nan(ntpts,2,nsims_sw);
    Ptot_sw_hl = nan(ntpts,nsims_sw);
    
    for jj_sim = 1:nsims_sw
        disp(['sample run ' num2str(jj_sim) ' of ' num2str(nsims_sw)]);
        
        U_sw_lh(:,:,jj_sim) = fi.sde([ssR(1),ssP(1)],tarr,muRsig_fit,muPsig_fit, ...
            ssP(3));
        Ptot_sw_lh(:,jj_sim) = U_sw_lh(:,2,jj_sim);
        
        U_sw_hl(:,:,jj_sim) = fi.sde([ssR(3),ssP(3)],tarr,muRsig_fit,muPsig_fit, ...
            ssP(1));
        Ptot_sw_hl(:,jj_sim) = U_sw_hl(:,2,jj_sim);
    end
    
    
    %% FIRST PASSAGE TIMES (PROTEIN VALUE)
    vals = linspace(0,ssP(3),25);
    
    tsw_lh = nan(nsims_sw,numel(vals));
    tsw_hl = nan(nsims_sw,numel(vals));
    
    for ii_val = 1:numel(vals)
        cur_val = vals(ii_val);
        for jj_sim = 1:nsims_sw
            cur_tsw_lh = tarr(find(Ptot_sw_lh(:,jj_sim) >= cur_val,1,'first'));
            if ~isempty(cur_tsw_lh)
                tsw_lh(jj_sim,ii_val) = cur_tsw_lh;
            end
            
            cur_tsw_hl = tarr(find(Ptot_sw_hl(:,jj_sim) <= cur_val,1,'first'));
            if ~isempty(cur_tsw_hl)
                tsw_hl(jj_sim,ii_val) = cur_tsw_hl;
            end
        end
    end
    
    
    figure
    errorbar(vals,nanmean(tsw_lh,1),nanstd(tsw_lh,[],1)/sqrt(nsims_sw),'ro');
    hold on;
    errorbar(vals,nanmean(tsw_hl,1),nanstd(tsw_hl,[],1)/sqrt(nsims_sw),'ko');
    hold off;
    xlabel('Ftz protein (au)')
    ylabel('first-passage time (min)')
    ylim([0,tmax])
    
    legend({'start at low state','start at high state'},'location','eastoutside')
    
    set_figure_defaults(gcf,'font');
    
    
    %% FIRST PASSAGE TIMES (RADIAL DISTANCE FROM STEADY STATE)
    vals = linspace(0,2e6,25);
    
    tsw_lh = nan(nsims_sw,numel(vals));
    tsw_hl = nan(nsims_sw,numel(vals));
    
    traj_distance_lh = sqrt( (squeeze(U_sw_lh(:,1,:)) - ssR(3)).^2 + ...
        (squeeze(U_sw_lh(:,2,:) - ssP(3)).^2) );
    traj_distance_hl = sqrt( (squeeze(U_sw_hl(:,1,:)) - ssR(1)).^2 + ...
        (squeeze(U_sw_hl(:,2,:) - ssP(1)).^2) );
    
    for ii_val = 1:numel(vals)
        cur_val = vals(ii_val);
        for jj_sim = 1:nsims_sw
            cur_tsw_lh = tarr(find(traj_distance_lh(:,jj_sim) <= cur_val,1,'first'));
            if ~isempty(cur_tsw_lh)
                tsw_lh(jj_sim,ii_val) = cur_tsw_lh;
            end
            
            cur_tsw_hl = tarr(find(traj_distance_hl(:,jj_sim) <= cur_val,1,'first'));
            if ~isempty(cur_tsw_hl)
                tsw_hl(jj_sim,ii_val) = cur_tsw_hl;
            end
        end
    end
    
    
    figure
    errorbar(vals,nanmean(tsw_lh,1),nanstd(tsw_lh,[],1)/sqrt(nsims_sw),'ro');
    hold on;
    errorbar(vals,nanmean(tsw_hl,1),nanstd(tsw_hl,[],1)/sqrt(nsims_sw),'ko');
    plot(sqrt((ssP(3)-ssP(1))^2 + (ssR(3)-ssR(1))^2)*[1 1],[0,tmax], ...
        '--','color',0.7*[1 1 1],'linewidth',1);
    hold off;
    xlabel('Euclidean distance (au)')
    ylabel('first-passage time (min)')
    ylim([0,tmax])
    
    legend({'start at low state','start at high state','dist btwn s.s.'}, ...
        'location','eastoutside')
    
    set_figure_defaults(gcf,'font');
end