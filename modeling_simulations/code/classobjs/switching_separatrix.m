classdef switching_separatrix < handle
    
    properties (SetAccess = private)
        filename;
        fi;
        
        Xq;
        Yq;
        Zq;
        
        predicted_on;
        false_positives;
        false_negatives;
        correctly_identified;
        frac_correctly_identified;
    end
    
    properties (Access = private)
        Mmax_pool;
        Mmin_pool;
    end
    
    methods
        
        % constructor
        function this = switching_separatrix(filename)
            if nargin > 0
                this.initialize(filename);
            end
        end
        
        function this = initialize(this,filename)
            this.filename = filename;
            ss = load(filename);
            
            this.Mmax_pool = ss.Mmax_pool;
            this.Mmin_pool = ss.Mmin_pool;
            this.fi = ss.fi;
            
            Xmax = this.Mmax_pool(:,1);
            Ymax = this.Mmax_pool(:,2);
            Zmax = this.Mmax_pool(:,3);
            
            Xmin = this.Mmin_pool(:,1);
            Ymin = this.Mmin_pool(:,2);
            Zmin = this.Mmin_pool(:,3);
            
            Xq1 = linspace(min([Xmax;Xmin]),max([Xmax;Xmin]),100);
            Yq1 = linspace(min([Ymax;Ymin]),max([Ymax;Ymin]),100);
            [this.Xq,this.Yq] = meshgrid(Xq1,Yq1);
            
            Zqmax = griddata(Xmax,Ymax,Zmax,this.Xq,this.Yq);
            Zqmin = griddata(Xmin,Ymin,Zmin,this.Xq,this.Yq);
            
            this.Zq = (Zqmin + Zqmax)/2;
        end
        
        function h = surf(this)
            h = surf(this.Xq,this.Yq,this.Zq);
            shading interp;
            
            set(h,'facealpha',0.5,'facecolor','k');
            xlabel('a_0'); ylabel('r_0'); zlabel('p_0');
            ylim([0,3e6]);
            view(-140,15);
        end
        
        % Determines whether given coordinates (x,y,z) are above the
        % separatrix.  Returns a boolean vector with one entry per
        % coordinate, true if that coordinate is above the surface.
        function [ix_above_surface,N_nearest,min_dist,min_ix] = ...
                above_surface(this,x,y,z)
            x = x(:);
            y = y(:);
            z = z(:);
            
            % get normal vectors to surface
            [Nx,Ny,Nz] = this.surfnorm();
            
            [min_dist,min_ix] = this.dist_from_surf(x,y,z);
            
            N_nearest = [Nx(min_ix)',Ny(min_ix)',Nz(min_ix)'];
            
            % dot product
            ix_above_surface = (diag(N_nearest*([x,y,z] - ...
                [this.Xq(min_ix)',this.Yq(min_ix)',this.Zq(min_ix)'])') > 0)';  
        end
        
        % Get normal vectors to surface.
        function [Nx,Ny,Nz] = surfnorm(this)
            [Nx,Ny,Nz] = surfnorm(this.Xq,this.Yq,this.Zq);
        end
        
        
        function [min_dist,min_ix] = dist_from_surf(this,x,y,z)
            x = x(:);
            y = y(:);
            z = z(:);
            
            all_dist = sqrt((x' - this.Xq(:)).^2 + (y' - this.Yq(:)).^2 + ...
                (z' - this.Zq(:)).^2);
            [min_dist,min_ix] = min(all_dist,[],1);
        end
        
        function [X,Y,Z] = displace_surface(this,dw,subsamp,zthresh,delete_nan)
            if nargin < 5
                delete_nan = false;
                if nargin < 4
                    zthresh = 0;
                    if nargin < 3
                        subsamp = 1;
                    end
                end
            end
            
            [Nx,Ny,Nz] = this.surfnorm();
            X = this.Xq + Nx*dw;
            Y = this.Yq + Ny*dw;
            Z = this.Zq + Nz*dw;
            
            X(X < 0) = NaN;
            Y(Y < 0) = NaN;
            Z(Z < zthresh) = NaN;
            
            if delete_nan
                ix_to_del = isnan(X) | isnan(Y) | isnan(Z);
                X(ix_to_del) = [];
                Y(ix_to_del) = [];
                Z(ix_to_del) = [];

                sample_ixs = randperm(numel(X),subsamp);
                X = X(sample_ixs);
                Y = Y(sample_ixs);
                Z = Z(sample_ixs);
%                 subsamp_freq = round(numel(X)/subsamp);
%                 X = X(1:subsamp_freq:end);
%                 Y = Y(1:subsamp_freq:end);
%                 Z = Z(1:subsamp_freq:end);
            else
                X = X(1:subsamp:end,1:subsamp:end);
                Y = Y(1:subsamp:end,1:subsamp:end);
                Z = Z(1:subsamp:end,1:subsamp:end);
            end
        end
        
        
        function [predicted_on,false_positives,false_negatives, ...
                correctly_identified,frac_correctly_identified] = ...
                prediction_accuracy(this,x,y,z,ground_truth_on,valid_ix,embryoid)
            if nargin < 7
                embryoid = ones(size(valid_ix));
            end
            
            predicted_on = this.above_surface(x,y,z);
            
            false_positives = (predicted_on & ~ground_truth_on);
            false_negatives = (~predicted_on & ground_truth_on);
            correctly_identified = (predicted_on & ground_truth_on) | ...
                (~predicted_on & ~ground_truth_on);
            frac_correctly_identified = sum(correctly_identified)/numel(ground_truth_on);
            
            disp(['prediction matches ' num2str(sum(correctly_identified)) ...
                ' out of ' num2str(sum(valid_ix)) ' nuclei (' ...
                num2str(100*frac_correctly_identified) '%) from ' ...
                num2str(numel(unique(embryoid(valid_ix)))) ' embryos']);
            disp(['- false positives: ' num2str(sum(false_positives)) ' (' ...
                num2str(100*sum(false_positives)/numel(ground_truth_on)) '%)']);
            disp(['- false negatives: ' num2str(sum(false_negatives)) ' (' ...
                num2str(100*sum(false_negatives)/numel(ground_truth_on)) '%)']);
            nembryos = numel(unique(embryoid));
            if nembryos > 1
                for ii = 1:nembryos
                    cie = sum(correctly_identified(embryoid(valid_ix) == ii));
                    nne = sum(valid_ix(embryoid == ii));
                    disp(['  embryo ' num2str(ii) ': prediction matches ' num2str(cie) ...
                        ' out of ' num2str(nne) ' nuclei (' num2str(100*cie/nne) '%)']);
                end
            end
            
            this.predicted_on = predicted_on;
            this.false_positives = false_positives;
            this.false_negatives = false_negatives;
            this.correctly_identified = correctly_identified;
            this.frac_correctly_identified = frac_correctly_identified;
        end
        
        
    end
    
end