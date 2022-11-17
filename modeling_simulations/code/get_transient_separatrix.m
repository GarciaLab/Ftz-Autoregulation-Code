% see Sootla et al. (2015) "Shaping pulses to control bistable systems"

clc; clear; close all;

tmax = 20;              % max simulation time (= time before gastrulation that autoregulation begins)
Nmax = 100;             % max # iterations for bisection

input_a0_range = [1e2,6e5];   % range of initial mRNA production rates
input_r0_range = [1e3,5e6];     % range of initial mRNA concentrations
input_p0_range = [0,1e7];       % range of initial protein concentrations

fi = ftz_instance();

[ssR,ssP] = fi.get_steady_states(inf);
if numel(ssR) ~= 3
    error('system is not bistable at infinite time');
end

prot_on_thresh = fi.prot_on_thresh;


%% ----- INITIAL EARLY mRNA PRODUCTION RATE + LEVEL, PROTEIN LEVEL ----- %%
Ne = 4;     % number random samples per dimension
Na0 = 24;	% number samples of a0

tol = 1e-3; % limit of error between min protein that switches and max that
            % does not switch (bisection)

a0_set = linspace(input_a0_range(1),input_a0_range(2),Na0);
Mmax_pool = [];
Mmin_pool = [];

for ii_a0 = 1:Na0
    disp(['Now processing ' num2str(ii_a0) ' of ' num2str(Na0) '...']);
    
    fi.set('a0',a0_set(ii_a0));
    
    r0_set = input_r0_range;
    p0_set = nan(size(r0_set));
    
    tic
    % seed with bisection
    for ii_r0 = 1:numel(r0_set)
        cur_input_p0 = input_p0_range(end);
        
        est_min_that_switches = input_p0_range(end);
        est_max_that_does_not_switch = input_p0_range(1);
        
        for ii_bis = 1:Nmax
            % set new input
            fi.set('r0',r0_set(ii_r0),'p0',cur_input_p0);
            
            % simulate
            [M,P,t] = fi.simulate([0,0],[0,tmax],'no_plot','no_save','off');
            Ptot = P(end) + fi.p(tmax);
            
            if isnan(M(end))
                error('simulation returned NaN');
            end
            
            % check if switch succeeded
            if Ptot > prot_on_thresh
                est_min_that_switches = cur_input_p0;
            elseif ii_bis == 1
                disp('maximum p0 is too low to switch system for minimum r0');
                break;
            else
                est_max_that_does_not_switch = cur_input_p0;
            end
            
            % update and continue
            cur_input_p0 = est_max_that_does_not_switch + ...
                (est_min_that_switches - est_max_that_does_not_switch)/2;
            
            if est_min_that_switches - est_max_that_does_not_switch < tol
                break
            end
        end
        p0_set(ii_r0) = cur_input_p0;
    end
    
    % initialize sets
    Mmax = [r0_set([1,numel(r0_set)])',p0_set([1,numel(p0_set)])'];
    Mmin = Mmax;
    
    cur_r0_range = r0_set([1,numel(r0_set)])';
    cur_p0_range = p0_set([1,numel(p0_set)])';
    
    N = 800;
    Npar = 2*Ne;
    Ntot = round(N/Npar);
    
    for ii = 1:Ntot
        % generate 2*Ne random samples
        rr0 = unifrnd(min(cur_r0_range),max(cur_r0_range),Ne,1);
        
        [~,iamax,~] = unique(Mmax(:,1));
        [~,iamin,~] = unique(Mmin(:,1));
        
        rr0_maxp0 = interp1(Mmax(iamax,1),Mmax(iamax,2),rr0,'prev');
        rr0_minp0 = interp1(Mmin(iamin,1),Mmin(iamin,2),rr0,'next');
        rr0_p0 = unifrnd(rr0_minp0,rr0_maxp0,Ne,1);
        
        if abs(diff(cur_p0_range)) > 0
            rp0 = unifrnd(min(cur_p0_range),max(cur_p0_range),Ne,1);
            
            [~,iamax,~] = unique(Mmax(:,2));
            [~,iamin,~] = unique(Mmin(:,2));
        
            rp0_maxr0 = interp1(Mmax(iamax,2),Mmax(iamax,1),rp0,'prev');
            rp0_minr0 = interp1(Mmin(iamin,2),Mmin(iamin,1),rp0,'next');
            rp0_r0 = unifrnd(rp0_minr0,rp0_maxr0,Ne,1);
        else
            rp0_r0 = [];
            rp0 = [];
        end
        
        % new inputs
        nin = [rr0,rr0_p0; ...
            rp0_r0,rp0];
        
        for jj = 1:size(nin,1)
            if ~isnan(nin(jj,2)) && ~isnan(nin(jj,1))
                % set new input
                fi.set('r0',nin(jj,1),'p0',nin(jj,2));
                
                % simulate
                [M,P,t] = fi.simulate([0,0],[0,tmax],'no_plot','no_save','off');
                Ptot = P(end) + fi.p(tmax);
                
                if isnan(M(end))
                    error('simulation returned NaN');
                end
                
                % check if switch succeeded
                if Ptot > prot_on_thresh
                    Mmax = [Mmax;nin(jj,:)];
                    Mmax = sortrows(Mmax,1,'ascend');   % above
                    
                    % prune
                    max_to_del = false(size(Mmax,1),1);
                    for kk = 1:size(Mmax,1)
                        if any(Mmax(1:kk-1,2) < Mmax(kk,2))
                            max_to_del(kk) = true;
                        end
                    end
                    Mmax = Mmax(~max_to_del,:);
                else
                    Mmin = [Mmin;nin(jj,:)];
                    
                    Mmin = sortrows(Mmin,1,'ascend');
                    min_to_del = false(size(Mmin,1),1);
                    for kk = 1:size(Mmin,1)
                        if any(Mmin(kk+1:end,2) > Mmin(kk,2))
                            min_to_del(kk) = true;
                        end
                    end
                    Mmin = Mmin(~min_to_del,:);
                end
            end
        end
    end
    toc
    
    Mmax_pool = [Mmax_pool;ones(size(Mmax,1),1)*a0_set(ii_a0),Mmax];
    Mmin_pool = [Mmin_pool;ones(size(Mmin,1),1)*a0_set(ii_a0),Mmin];
end

figure
scatter3(Mmax_pool(:,1),Mmax_pool(:,2),Mmax_pool(:,3)); hold on;
scatter3(Mmin_pool(:,1),Mmin_pool(:,2),Mmin_pool(:,3)); hold off;
xlabel('a0'); ylabel('r0'); zlabel('p0');


%% ----- SAVE RESULT ----- %%
save('trans_switching_separatrix_a0-r0-p0','Mmax_pool','Mmin_pool','fi');