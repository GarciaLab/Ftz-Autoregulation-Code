function set_figure_defaults(fig,varargin)
set_font = false;
set_lines = false;
linewidth = 2;

if isempty(varargin)
    varargin = {'both'};
end

v = 1;
while v <= numel(varargin)
    switch varargin{v}
        case 'font'
            set_font = true;
        case 'lines'
            set_lines = true;
        case 'both'
            set_font = true;
            set_lines = true;
        case 'linewidth'
            assert(numel(varargin) > v);
            linewidth = varargin{v+1};
            v = v+1;
        otherwise
            error('unsupported specification %s',varargin{v});
    end
    v = v + 1;
end

ax = fig.Children;

for ii = 1:length(ax)
    if set_font
        set(ax(ii),'FontSize',14);
        set(ax(ii),'LineWidth',1);
        set(ax(ii),'box','off');
    end
    
    if set_lines
        h = ax(ii).Children;
        for jj = 1:length(h)
            if ~strcmp(class(h(jj)),'matlab.graphics.primitive.Image')
                set(h(jj),'LineWidth',linewidth);
            end
        end
    end
end

end