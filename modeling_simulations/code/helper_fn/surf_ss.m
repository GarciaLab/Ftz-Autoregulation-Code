function [h,ss_set_updated] = surf_ss(ss_set,color_set)
if nargin < 2 || isempty(color_set)
    color_set = {'k','b','g','r'};
end

for ii = 1:numel(ss_set)
    if iscell(ss_set)
        ss = ss_set{ii};
    elseif isstruct(ss_set)
        ss = ss_set;
    else
        error('separatrix must be a struct or cell of structs');
    end
    
    Xmax = ss.Mmax_pool(:,1);
    Ymax = ss.Mmax_pool(:,2);
    Zmax = ss.Mmax_pool(:,3);
    
    Xmin = ss.Mmin_pool(:,1);
    Ymin = ss.Mmin_pool(:,2);
    Zmin = ss.Mmin_pool(:,3);
    
    Xq1 = linspace(min([Xmax;Xmin]),max([Xmax;Xmin]),100);
    Yq1 = linspace(min([Ymax;Ymin]),max([Ymax;Ymin]),100);
    [Xqt,Yqt] = meshgrid(Xq1,Yq1);
    
    Vqmax = griddata(Xmax,Ymax,Zmax,Xqt,Yqt);
    Vqmin = griddata(Xmin,Ymin,Zmin,Xqt,Yqt);
    
    Vqt = (Vqmin + Vqmax)/2;
    
    if ii == 1
        h = surf(Xqt,Yqt,Vqt);
    else
        h = [h,surf(Xqt,Yqt,Vqt)];
    end
    shading interp;
    hold on;
    
    ss.Xq = Xqt;
    ss.Yq = Yqt;
    ss.Vq = Vqt;
    
    if iscell(ss_set)
        ss_set_updated{ii} = ss;
    elseif isstruct(ss_set)
        ss_set_updated = ss;
    end
end
hold off;

for ii = 1:numel(ss_set)
    set(h(ii),'facealpha',0.5,'facecolor',color_set{ii});
end

xlabel('a_0'); ylabel('r_0'); zlabel('p_0');
ylim([0,3e6]);

end