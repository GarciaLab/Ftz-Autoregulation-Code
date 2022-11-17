% creator: Mindy Perkins

classdef ftz_instance < handle
    
    properties
        late_mRNA;
        protein;
        
        markersz = 50;
    end
    
    properties (SetAccess = private)
        input_type = 'early_protein';	% specifies whether input is early mRNA or early protein
        
        c = 0.45;               % i/o fn ratio endogenous/transgene
        alpha = 0.08165;        % translation rate (prot AU min^-1 (mRNA AU)^-1)
        tauR = 7;               % Ftz mRNA half-life A (min), our data, Edgar et al. (1986b)
        tauP = 7.9;             % Ftz protein half-life (min), Bothma et al. (2018)
        gammaR;                 % Ftz mRNA decay rate (mRNA AU min^-1)
        gammaP;                 % Ftz protein decay rate (prot AU min^-1)
        
        % early mRNA
        a0 = 0;                 % early mRNA production rate
        r0 = 0;                 % initial early mRNA
        p0 = 0;                 % initial early protein
        b = 0.0479;             % decay rate of early mRNA production
        r;                      % analytical model early mRNA
        p;                      % analytical model early protein
        c1;                     % for default analytical model early protein
        c2;
        c3;
        
        % input/output function parameters
        K = 1.216e6;
        a = 5.555e5;
        n = 3.264;
        f;
        l = 0;                  % leakiness
        
        time_pts;               % array of time_pts for trace measurements
        
        sim_tt;                 % simulated times points
        sim_late_mRNA;          % simulated late mRNA concentration
        sim_protein;            % simulated protein concentration
        
        sim_sensitivities;      % sensitivity analysis results for simulation
        
        force_bistable = false;
        prot_lower_bound = [];
        
        prot_on_thresh = 0.5e6;    % protein threshold for "on" fate
    end
    
    methods
        
        % Constructor.
        function this = ftz_instance()
            this.gammaP = log(2)/this.tauP;
            this.gammaR = log(2)/this.tauR;
            
            this.r = @(t) (this.a0/(this.gammaR - this.b))*exp(-this.b*t) + ...
                (this.r0 - this.a0/(this.gammaR - this.b))*exp(-this.gammaR*t);
            
            this.reset_p();
            
            this.f = @(t,x) (this.a*x.^this.n)./(this.K.^this.n + x.^this.n)+ this.l;
        end
        
        
        % Sets p to the default linear system model for input.
        function reset_p(this)
            this.set_analytical_p_params();
            
            this.p = @(t) this.c1*(this.alpha/((this.gammaR - this.b)*(this.gammaP-this.b)))*exp(-this.b*t) + ...
                this.c2*(this.alpha/(this.gammaP - this.gammaR))*exp(-this.gammaR*t) + ...
                this.c3*exp(-this.gammaP*t);
        end
        
        
        % Set object properties.
        function set(this,varargin)
            % parse optional inputs
            v = 1;
            while v <= numel(varargin)
                switch varargin{v}
                    case 'f'    % i/o fn
                        assert(numel(varargin) >= v);
                        v = v+1;
                        fp = varargin{v};
                        
                        if ~isa(fp,'function_handle')
                            error('f must be a function handle');
                        else
                            this.f = fp;
                        end
                    case 'r'    % early mRNA decay
                        assert(numel(varargin) >= v);
                        v = v+1;
                        rp = varargin{v};
                        
                        if ~isa(rp,'function_handle')
                            error('r must be a function handle');
                        else
                            this.r = rp;
                        end
                    case 'p'    % early protein input
                        assert(numel(varargin) >= v);
                        v = v+1;
                        pp = varargin{v};
                        
                        if ~isa(pp,'function_handle')
                            error('p must be a function handle');
                        else
                            if this.force_bistable
                                this.p = @(t) max(this.prot_lower_bound,pp(t));
                            else
                                this.p = pp;
                            end
                        end
                    case 'tauR' % mRNA half-life
                        assert(numel(varargin) >= v);
                        v = v+1;
                        this.tauR = varargin{v};
                        this.gammaR = log(2)/this.tauR;   % Ftz mRNA decay rate (AU min^-1)
                    case 'tauP' % protein half-life
                        assert(numel(varargin) >= v);
                        v = v+1;
                        this.tauP = varargin{v};
                        this.gammaP = log(2)/this.tauP;   % Ftz protein decay rate (AU min^-1)
                    case 'alpha' % translation rate
                        assert(numel(varargin) > v);
                        this.(varargin{v}) = varargin{v+1};
                        v = v+1;
                    otherwise
                        assert(numel(varargin) > v);
                        this.(varargin{v}) = varargin{v+1};
                        v = v+1;
                end
                v = v+1;
            end
            
            this.reset_simulated_results();
            this.set_analytical_p_params();
        end
        
        
        % Set the lump terms used in the analytical solution to the linear
        % system approximation of early mRNA and protein.  Irrelevant if
        % this.p has been manually set by user.
        function set_analytical_p_params(this)
            this.c1 = this.a0;
            this.c2 = this.r0 - (1/(this.gammaR - this.b))*this.c1;
            this.c3 = this.p0 - (this.alpha/(this.gammaP - this.gammaR))*this.c2 - ...
                (this.alpha/((this.gammaR - this.b)*(this.gammaP-this.b)))*this.c1;
        end
        
        
        % Sets flag to lower bound all uses of the input by the minimum
        % input value that guarantees bistability.
        function set_force_bistable(this)
            plb = this.get_input_monostable();
            
            if isnan(plb)
                error('no input values guarantee bistability');
            end
            
            this.force_bistable = true;
            this.prot_lower_bound = ceil(plb*1e-4)*10000;
            pp = this.p;
            this.set('p',pp);
        end
       
        
        % Clear any stored simulated results.
        function reset_simulated_results(this)
            this.sim_late_mRNA = [];
            this.sim_protein = [];
        end
        
        
        % Plot the input-output function for this system.
        function fig = plot_io(this,protein,tau_set)
            if nargin < 3
                tau_set = 0;	% times to plot (min)
                if nargin < 2
                    error('usage is: obj.plot_io(x,tau_set)');
                end
            end
            
            fig = figure;
            for ii_tau = 1:length(tau_set)
                h(ii_tau) = plot(protein, ...
                    this.c*this.f(tau_set(ii_tau),protein));
                hold on;
            end
            hold off;
            
            this.default_legend(h,tau_set);
            title({'Input/output function','(endogenous autoregulatory)'})
            xlabel('Ftz protein (AU)');
            ylabel('transcription rate (AU)');
            
            set_figure_defaults(fig);
        end
        
        
        % Numerically identify the steady states of the system by
        % evaluating crossings of the transfer function with a line of
        % slope 1 at time tau.  Optionally replace specified variables with
        % specified values.
        function [ssR,ssP,T] = get_steady_states(this,tau,varargin)
            c = this.c;
            gammaR = this.gammaR;
            alpha = this.alpha;
            gammaP = this.gammaP;
            
            calc_a1 = true;
            calc_a2 = true;
            
            % parse optional inputs
            v = 1;
            while v <= numel(varargin)
                switch varargin{v}
                    case 'c'
                        assert(numel(varargin) > v);
                        v = v+1;
                        c = varargin{v};
                    case 'alpha'
                        assert(numel(varargin) > v);
                        v = v+1;
                        alpha = varargin{v};
                    case 'gammaR'
                        assert(numel(varargin) > v);
                        v = v+1;
                        gammaR = varargin{v};
                    case 'gammaP'
                        assert(numel(varargin) > v);
                        v = v+1;
                        gammaP = varargin{v};
                    case 'tauR'
                        assert(numel(varargin) > v);
                        v = v+1;
                        gammaR = log(2)/varargin{v};
                    case 'tauP'
                        assert(numel(varargin) > v);
                        v = v+1;
                        gammaP = log(2)/varargin{v};
                    case 'a1'
                        calc_a1 = false;
                        assert(numel(varargin) > v);
                        v = v+1;
                        a1 = varargin{v};
                    case 'a2'
                        calc_a2 = false;
                        assert(numel(varargin) > v);
                        v = v+1;
                        a2 = varargin{v};
                    otherwise
                        error('unsupported option %s',varargin{v});
                end
                v = v+1;
            end
            
            if numel(tau) > 1
                error('tau must be a scalar');
            end
            
            if calc_a1
                a1 = c/gammaR;
            end
            
            if calc_a2
                a2 = alpha/gammaP;
            end
            
            T = @(R) a1* ...
                    this.f(tau,a2*R + this.p(tau));

            TT = @(R) T(R(:)) - R(:);
            
            Rlims = [0,T(1e20) + 1]; % search until saturation of TT
            
            rr = linspace(Rlims(1),Rlims(2),1000);
            ixs = [(diff(sign(TT(rr))) ~= 0);false];
            guesses = rr(ixs);
            
            ssR = find_bistable_xings(TT,guesses);
            ssP = a2*ssR;
        end
        
        
        % Numerically calculate the minimum protein input value for which
        % the system is bistable.
        function cur_p = get_input_monostable(this,pmin,pmax,tol,Nmax)
            if nargin < 5
                Nmax = 100;
                if nargin < 4
                    tol = 1e-2;
                    if nargin < 3
                        pmax = 1e5;
                        if nargin < 2
                            pmin = 0;
                        end
                    end
                end
            end
            
            this.set('p',@(t) pmin);
            [ssR,~] = this.get_steady_states(0);
            if numel(ssR) > 1
                cur_p = pmin;
                return;
            end
            
            % first lower pmax until enter a bistable region
            while true
                this.set('p',@(t) pmax);
                [ssR,~] = this.get_steady_states(0);
                if numel(ssR > 1)
                    break;
                end
                
                pmax = 0.9*pmax;
            end
            
            % now use bisection
            cur_p = pmax;
            for ii = 1:Nmax
                this.set('p',@(t) cur_p);
                [ssR,~] = this.get_steady_states(0);
                
                if numel(ssR) == 3
                    pmax = cur_p;
                else
                    pmin = cur_p;
                end
                cur_p = pmin + (pmax - pmin)/2;
                
                if (pmax - pmin) < tol
                    break
                end
                
                if ii == Nmax
                    disp('Reached max number iterations--terminating...')
                end
            end
            
            this.reset_p();
        end
            
        
        % Plot a graphical test to determine bistability from
        % intersections of the input/output function with a line (see
        % write-up).
        function fig = plot_ss_graphical_test(this,xx,tau_set,opt,plotopt)
            plot_protein = true;
            new_fig = false;
            
            if nargin < 5
                plotopt = 'no_new_fig';
                if nargin < 4
                    opt = 'protein';
                    if nargin < 3
                        error('usage is: obj.plot_ss_graphical_test(RR,tau_set)');
                    end
                end
            end
            
            if strcmp(opt,'rna')
                plot_protein = false;
            end
            
            if strcmp(plotopt,'new_fig')
                new_fig = true;
            end
            
            ntau = length(tau_set);
            cl = lines(ntau);
            
            if new_fig
                fig = figure;
            end
            
            if plot_protein
                h = plot(xx,this.gammaP*xx,'--','Color',0.7*[1 1 1]);
            else
                h = plot(xx,this.gammaR*xx,'--','Color',0.7*[1 1 1]);
            end
            hold on;
            for ii_tau = 1:ntau
                [ssR,ssP,T] = this.get_steady_states(tau_set(ii_tau));
                
                ssT = T(ssR);
                
                if plot_protein
                    h(ii_tau+1) = plot(xx, ...
                        this.alpha*T((this.gammaP/this.alpha)*xx),'-','Color', ...
                        cl(ii_tau,:));
                else
                    h(ii_tau+1) = plot(xx,this.gammaR*T(xx),'-','Color', ...
                        cl(ii_tau,:));
                end
                
                % check stability of identified points
                fprime = @(x) ((this.a*this.n*x.^(this.n-1))*this.K^this.n)./((this.K^this.n + x.^this.n).^2);
                lambda_pos = (-(this.gammaP + this.gammaR) + sqrt((this.gammaP + this.gammaR)^2 - ...
                    4*(this.gammaR*this.gammaP - this.alpha*this.c*fprime(ssP))))/2;
                lambda_neg = (-(this.gammaP + this.gammaR) - sqrt((this.gammaP + this.gammaR)^2 - ...
                    4*(this.gammaR*this.gammaP - this.alpha*this.c*fprime(ssP))))/2;
                
                s_ix = (lambda_pos < 0) & (lambda_neg < 0); % ix of stable pts
                us_ix = ~s_ix;  % ix of unstable pts
                
                if plot_protein
                    scatter(ssP(s_ix),this.alpha*ssT(s_ix), ...
                        this.markersz,'o','LineWidth',1, ...
                        'MarkerFaceColor',cl(ii_tau,:), ...
                        'MarkerEdgeColor',cl(ii_tau,:));
                    scatter(ssP(us_ix),this.alpha*ssT(us_ix), ...
                        this.markersz,'o','LineWidth',1, ...
                        'MarkerEdgeColor',cl(ii_tau,:));
                else
                    scatter(ssR(s_ix),this.gammaR*ssT(s_ix), ...
                        this.markersz,'o','LineWidth',1, ...
                        'MarkerFaceColor',cl(ii_tau,:), ...
                        'MarkerEdgeColor',cl(ii_tau,:));
                    scatter(ssR(us_ix),this.gammaR*ssT(us_ix), ...
                        this.markersz,'o','LineWidth',1, ...
                        'MarkerEdgeColor',cl(ii_tau,:));
                end
            end
            hold off;
            
            % generate legend
            if plot_protein
                this.default_legend(h,tau_set, ...
                    'prior_entries',{'$$\gamma_P P_\tau^*$$'});
                
                title('$$\frac{\alpha c}{\gamma_R}f\left(\tau,P_\tau^* + p(\tau)\right)$$','interpreter','latex');
                xlabel('protein $$P_\tau^*$$','interpreter','latex');
            else
                this.default_legend(h,tau_set, ...
                    'prior_entries',{'$$\gamma_R R_\tau^*$$'});
                
                title('$$cf\left(\tau,\frac{\alpha}{\gamma_P}R_\tau^* + p(\tau)\right)$$','interpreter','latex');
                xlabel('mRNA concentration $$R_\tau^*$$','interpreter','latex');
            end
            
            set_figure_defaults(gcf);
        end
        
        
        % Simulate the dynamical system.  Option to calculate the
        % parameter sensitivities df/dtheta following the method in
        % Dickinson and Gelinas (2007).
        %
        % Inputs:
        % - ic, 2-element vector of initial conditions [mRNA,protein]
        % - tt, time points for which to return simulation results
        % - plot_opt, 'plot' if function should show plot
        % - save_opt, 'save' if results should be stored in this object
        % - sens_opt, 'on' if parameter sensitivities should be calculated
        function [M,P,t,Z] = simulate(this,ic,tt,plot_opt,save_opt,sens_opt)
            if nargin < 6
                sens_opt = 'off';
                if nargin < 5
                    save_opt = 'no_save';
                    if nargin < 4
                        plot_opt = 'no_plot';
                        if nargin < 3
                            tt = [];
                            if nargin < 2
                                ic = [];
                            end
                        end
                    end
                end
            end
            
            if strcmp(plot_opt,'plot')
                gen_plot = true;
            else
                gen_plot = false;
            end
            
            if strcmp(save_opt,'save')
                save_res = true;
            else
                save_res = false;
            end
            
            if strcmp(sens_opt,'on')
                sim_sens = true;
                disp('SENSITIVITY ANALYSIS CURRENTLY DOES NOT SUPPORT TIME-DEPENDENT I/O FN');
            else
                sim_sens = false;
            end
            
            if isempty(ic)
                ic = [this.late_mRNA(1),this.protein(1)];
            end
            
            if isempty(tt)
                tt = this.time_pts;
            end
            
            if ~isempty(this.sim_late_mRNA) && isequal(this.sim_tt,tt)
                disp(['returning precalculated simulated results ', ...
                    '(to overwrite first call this.reset_simulated_results)']);
                
                M = this.sim_late_mRNA;
                P = this.sim_protein;
            else
                if sim_sens
                    % number parameters for sensitivity (excludes i.c.s)
                    nparam = 7;
                    % order: c, a, K, n, gammaR, alpha, gammaP
                    
                    % for sensitivity analysis
                    dfdK = @(t,x) -this.f(t,x) .* this.n*this.K^(this.n-1) ./ ...
                        (this.K^this.n + x.^this.n);
                    dRdn = @(t,x) this.f(t,x) .* ...
                        (this.K^this.n).*log(x./this.K)./(this.K^this.n + x.^this.n);
                    dodedth = @(t,R,P) [ ...
                        this.f(t,P+this.p(t)); 0; ...
                        (this.c/this.a)*this.f(t,P+this.p(t)); 0; ...
                        this.c*dfdK(t,P+this.p(t)); 0; ...
                        this.c*dRdn(t,P+this.p(t)); 0; ...
                        -R; 0; 0; R; 0; -P];
                    dfdx = @(x) this.a*this.n*x.^(this.n-1)*this.K^this.n ./ ...
                        (this.K^this.n + x.^this.n).^2; %  DOES NOT SUPPORT TIME-DEPENDENT I/O
                    J = @(x) [-this.gammaR, this.c*dfdx(x); ...
                        this.alpha, -this.gammaP];
                    
                    odefun = @(t,X) ...
                        [this.c*this.f(t,X(2) + this.p(t)) - this.gammaR*X(1); ...
                        this.alpha*X(1) - this.gammaP*X(2); ...
                        dodedth(t,X(1),X(2)) + kron(eye(nparam),J(X(2)+this.p(t)))*X(3:end)];
                    ic = [ic,zeros(1,2*nparam)];
                else % w/o sensitivity calculations
                    odefun = @(t,X) ...
                        [this.c*this.f(t,X(2) + this.p(t)) - this.gammaR*X(1); ...
                        this.alpha*X(1) - this.gammaP*X(2)];
                end
                
                [t,X] = ode45(@(t,x) odefun(t,x),tt,ic);
                
                M = X(:,1);
                P = X(:,2);
                Z = X(:,3:end);
                
                if save_res
                    this.sim_late_mRNA = M;
                    this.sim_protein = P;
                    this.sim_tt = tt';
                    this.sim_sensitivities = Z;
                end
            end
            
            if gen_plot
                figure
                plot(M,P);
            end
        end
        
        
        % Stochastic simulation using Euler-Maruyama method.
        %
        % Inputs:
        % - uic, 2-element vector of initial conditions [late mRNA, total
        %   protein]
        % - tarr, time points for which to return simulation results
        % - muRsig_fit, function or fit object with standard deviation of
        %   mRNA noise as function of mRNA concentration
        % - muPsig_fit, function or fit object with standard deviation of
        %   protein noise as function of protein concentration
        % - xing_thresh, protein value that, if crossed, triggers
        %   simulation end
        function [u,tarr] = sde(this,uic,tarr,muRsig_fit,muPsig_fit, ...
                xing_thresh)
            if nargin < 6
                xing_thresh = [];
            end
            
            trigger_stop = false;
            greater_than = false;
            if ~isempty(xing_thresh)
                trigger_stop = true;
                if xing_thresh > uic(2)
                    greater_than = true;
                end
            end
            
            dt = tarr(2)-tarr(1);
            ntpts = numel(tarr);
            
            r_early = @(t) (1/(this.gammaR - this.b))*this.c1*exp(-this.b*t) + ...
                this.c2*exp(-this.gammaR*t);
            
            u = nan(2,ntpts);
            
            % Rlate, Ptot
            u(:,1) = uic';
            dW = sqrt(dt)*randn(2,ntpts-1);
            
            odefun = @(t,X) ...
                [this.c*this.f(t,max(X(2),this.p(inf))) - this.gammaR*X(1); ...
                this.alpha*(X(1)+r_early(t)) - this.gammaP*X(2)];
            
            % Euler-Maruyama method
            for ii = 2:ntpts
                u(:,ii) = u(:,ii-1) + odefun(tarr(ii),u(:,ii-1))*dt + ...
                    [1;1].*[muRsig_fit(u(1,ii-1))*dW(1,ii-1); ...
                    muPsig_fit(u(2,ii-1))*dW(2,ii-1)];
                
                % floor concentration at 0
                u(:,ii) = max(u(:,ii),0);
                
                if trigger_stop
                    if greater_than && (u(2,ii) > xing_thresh)
                        break;
                    elseif ~greater_than && (u(2,ii) <= xing_thresh)
                        break;
                    end
                end
            end
            
            u = u';
        end
        
        
        % Calculate the Hessian for square error integrated over time and
        % normalized by species maximum.  Implicitly treats the parameters
        % of this instance as the nominal.
        function [H,params] = sensitivity_hessian(this)
            if isempty(this.sim_sensitivities)
                error('cannot calculate: sim_sensitivities does not exist');
            end
            
            Ns = 2;  % number of species
            nparam = round(size(this.sim_sensitivities,2)/Ns);
            Tc = this.sim_tt(end) - this.sim_tt(1);  % time span
            
            % mRNA
            mRNA_max = max(this.sim_late_mRNA);
            
            % protein
            prot_max = max(this.sim_protein);
            
            temp = (this.sim_sensitivities.'*this.sim_sensitivities);
            
            H = (1/(Ns*Tc))*((1/mRNA_max^2)*temp(1:nparam,1:nparam) + ...
                (1/prot_max^2)*temp(nparam+1:end,nparam+1:end));
            
            % order: c, a, K, n, gammaR, alpha, gammaP
            params = [this.c;this.a;this.K;this.n;this.gammaR;this.alpha;this.gammaP];
        end
        
        
        % Generate a default plot legend in which each line is labeled by
        % the time and (optionally) early mRNA level to which it
        % corresponds.
        function l = default_legend(this,h,tau_set,varargin)
            prior_entries = {};
            
            % parse optional inputs
            v = 1;
            while v <= numel(varargin)
                switch varargin{v}
                    case 'prior_entries'
                        assert(numel(varargin) >= v);
                        v = v+1;
                        prior_entries = varargin{v};
                    otherwise
                        error('unsupported option %s',varargin{v});
                end
                v = v+1;
            end
            
            legend_entries = cell(numel(tau_set)+numel(prior_entries),1);
            legend_entries(1:numel(prior_entries)) = prior_entries;
            
            l2txt = @(t) ['$$\tau = ' num2str(t),'$$'];
            
            legend_entries(numel(prior_entries)+1:end) = ...
                cellfun(l2txt,num2cell(tau_set),'UniformOutput',false);
            l = legend(h,legend_entries,'interpreter','latex');
            set(legend,'location','southeastoutside','box','off');
        end
        
        
        % Scatterplot instantaneous steady-state mRNA vs. steady-state
        % protein.  Helper function for public plotting methods.
        function h = scatter_ss(this,tau,markercolor)
            if nargin < 3
                markercolor = 'k';
            end
            
            if length(markercolor) > 3
                markeralph = markercolor(4);
                markercolor = markercolor(1:3);
            else
                markeralph = 1;
                markercolor = markercolor(1:3);
            end
            
            [ssR,ssP,~] = this.get_steady_states(tau);
            if numel(ssR) == 3 % bistable
                h = scatter(ssR([1,3]),ssP([1,3]),this.markersz,'o', ...
                    'LineWidth',1,'MarkerFaceColor',markercolor, ...
                    'MarkerEdgeColor',markercolor, ...
                    'MarkerFaceAlpha',markeralph, ...
                    'MarkerEdgeAlpha',markeralph);
                hold on;
                scatter(ssR(2),ssP(2),this.markersz,'o', ...
                    'LineWidth',1,'MarkerEdgeColor',markercolor, ...
                    'MarkerFaceAlpha',markeralph, ...
                    'MarkerEdgeAlpha',markeralph);
                hold off;
            else % monostable
                h = scatter(ssR,ssP,this.markersz,'o','LineWidth',1, ...
                    'MarkerFaceColor',markercolor, ...
                    'MarkerEdgeColor',markercolor, ...
                    'MarkerFaceAlpha',markeralph, ...
                    'MarkerEdgeAlpha',markeralph);
            end
        end
        
        
    end
    
end