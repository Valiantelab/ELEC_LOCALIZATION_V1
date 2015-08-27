function [map, limits, cmap]=vals2Colormap(vals,type,cmap,minmax)
%function [map, limits]=vals2Colormap(vals,type,cmap,minmax)
% creates colormap (i.e., colorscale)
%
% Required Inputs:
%  vals - vector of values
%  type - 'absmax', 'minmax', 'justpos', 'justneg', or a 2D numeric vector.
%  If a 2D numeric vector, the first element is taken as the min and the
%  second as the max value of the colormap.
%  cmap - 'parula', 'jet' or 'rb'. No effect if 'justpos' or 'justneg' selected as type.
%
%  Note 'rb' is still a work in progress!
%
% Optional Input:
%  minmax - Two element vector specifying the low and high limits of the
%           color scale (in units of vals). {default: derived from vals}
%
% Author: David Groppe
% Mehtalab 2012
%
% Future work: get rid of minmax argument. I think it's only useful when
% type is justpos or justneg.

% History:
% 2015-4: Made compatible with Matlab 2014's parula cmap: DG
% 2015-4: Can handle saturated colormaps: Kathrin Müsch

if nargin<3,
    if verLessThan('matlab','8.0.1')
        cmap='jet';
    else
        cmap='parula';
    end
end

if nargin<4
    minmax=[];
end

if isnumeric(type)
   minmax=type;
   type='minmax';
end


if sum(isnan(vals))
    error('You have one or more NaN values in "vals" variable.');
elseif sum(isinf(vals))
    error('You have one or more inf values in "vals" variable.');
end

if isnumeric(type)
    %manually specified minmax scaling
    if length(type)~=2,
        error('Numeric value of "type" must be two dimensional.');
    end
    switch cmap
        case {'parula'}
            rgb_vals=colormap('parula');
        case {'jet'}
            rgb_vals=colormap('jet');
            %load jet_cmap
        case {'rb'}
            %% work in progress? ??
            rgb_vals=[linspace(0,1,128) ones(1,128); ...
                linspace(0,1,128) linspace(1,0,128); ...
                ones(1,128) linspace(1,0,128)]';
        otherwise
            error('I do not recognize cmap of type: %s',cmap);
    end
    n_colors=size(rgb_vals,1);
    
    cbar_max=max(type);
    cbar_min=min(type);
    temp_vals=vals-cbar_min;
    temp_vals=temp_vals/(cbar_max-cbar_min); %should range from 0 to 1
    temp_vals=temp_vals*(n_colors-1)+1; %should range from 1 to n_colors
    vals_col=round(temp_vals);
    map=rgb_vals(vals_col,:);

elseif strcmpi(type,'absmax')
    %absmax color scaling
    if verLessThan('matlab','8.0.1')
        rgb_vals=colormap('jet');
        %load jet_cmap. I did this in the past instead of calling colormap.
        %I don't remember why.
    else
        rgb_vals=colormap('parula');
    end
    n_colors=size(rgb_vals,1);
    
    if isempty(minmax)
        cbar_max=max(abs(vals));
        cbar_min=-cbar_max;
    else
        cbar_max=max(abs(minmax));
        cbar_min=-cbar_max;
    end
    temp_vals=((vals/cbar_max)+1)/2; %should range from 0 to 1
    temp_vals=round(temp_vals*(n_colors-1)+1); %should range from 1 to n_colors
    map=rgb_vals(temp_vals,:);
elseif strcmpi(type,'minmax')
    %minmax color scaling
    if verLessThan('matlab','8.0.1')
        rgb_vals=colormap('jet');
        %load jet_cmap
    else
        rgb_vals=colormap('parula');
    end
    n_colors=size(rgb_vals,1);
    
    if isempty(minmax)
        cbar_max=max(vals);
        cbar_min=min(vals);
    else
        cbar_max=max(minmax);
        cbar_min=min(minmax);
    end
    
    cval = linspace(cbar_min,cbar_max,n_colors)';
    nVals=length(vals);
    vals_col=zeros(nVals,1);
    for myidx = 1:nVals,
        vals_col(myidx) = nearest(cval,vals(myidx))';
    end
    map=rgb_vals(vals_col,:);
    
elseif strcmpi(type,'justpos') || strcmpi(type,'justneg')
    %colormap by limits
    if strcmpi(type,'justpos');
        % Just positive values
        %rgb_vals=colormap('autumn');
        load autumn_cmap
        cmap='autumn';
        if isempty(minmax)
            cbar_max=max(vals);
            cbar_min=0;
        else
            cbar_max=max(minmax);
            cbar_min=min(minmax);
        end
    else
        % Just positive values
        %rgb_vals=colormap('winter');
        load winter_cmap
        cmap='winter';
        if isempty(minmax)
            cbar_max=0;
            cbar_min=min(vals);
        else
            cbar_max=max(minmax);
            cbar_min=min(minmax);
        end
    end
    
    % Computer color indices for data values
    n_colors=size(rgb_vals,1);
    delt=(cbar_max-cbar_min)/(n_colors-1);
    data2color_values=cbar_min:delt:cbar_max;
    n_vals=length(vals);
    map=zeros(n_vals,3);
    for i=1:n_vals
        color_id=findTpt(vals(i),data2color_values);
        map(i,:)=rgb_vals(color_id,:);
    end
else
    error('Invalid value for ''type''.');
end

limits=[cbar_min cbar_max];