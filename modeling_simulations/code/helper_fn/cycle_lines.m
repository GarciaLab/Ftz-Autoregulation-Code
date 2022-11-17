function cycle_lines(h,varargin)

styles = false;
colors = true;

% parse optional inputs
v = 1;
while v <= numel(varargin)
    switch varargin{v}
        case 'styles'
            styles = true;
        case 'colors'
            colors = true;
        case 'nocolors'
            colors = false;
        case 'nostyles'
            styles = false;
        otherwise
            error('unsupported option %s',varargin{v});
    end
    v = v+1;
end

if colors
    c = num2cell(get(0,'DefaultAxesColorOrder'),2);
    set(h,{'Color'},c(rem((1:numel(h))-1,numel(c))+1));
end

if styles
    l = cellstr(get(0,'DefaultAxesLineStyleOrder'));
    set(0,'DefaultAxesLineStyleOrder','-|--|:|-.');
    set(h,{'LineStyle'},l(rem((1:numel(h))-1,numel(l))+1));
end

end