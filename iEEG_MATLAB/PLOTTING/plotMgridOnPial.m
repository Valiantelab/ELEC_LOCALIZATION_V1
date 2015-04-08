function plotMgridOnPial(fsub,hem,printEm,gridStems,gridDims)
%function plotMgridOnPial(fsub,hem,printEm,gridStems,gridDims)
% 
% Makes romni or lomni plots on a gray pial and DK atlas pial of a
% patient's brain with electrodes colored according to mgrid file colors.
%
% Inputs:
%  fsub  - FreeSurfer subject directory
%  hem   - ['r' or 'l'] Hemisphere to plot
%  printEm - If nonzero jpgs of the figures will be made in the subject's
%            elec_recon folder
%  gridStems - A string or cell array of strings indicating the stems of
%              grids (e.g., 'LGr' or {'LGrA','LGrB'}
%  gridDims - A n x 2 array where n=length(gridStems). Each row of the array 
%             indicates the number of rwos and columns in the grid (e.g.,
%             [8 8] or [4 8; 4 8]);
%
% Examples:
% % Strips-only patient
% plotMgridOnPial('TWH10','l',0);
%
% % Grid, strips, & depth patient
% plotMgridOnPial('TWH11','l',0,'LTGrid',[8 8]);
%
%
% Author:
% David M. Groppe
% March, 2015
%
% Future Work:
%  -Check to see if this really does work with grids that have been cut
%  -Add an option to make depths invisible
%  -Currently I use derive_grid_lines.m to find the neighbors of
%  electrodes. It might be easier to get this directly from the mgrid file.

if nargin<4
    gridStems=[];
elseif ~iscell(gridStems)
    %convert string to a singleton cell array
    temp=gridStems;
    clear gridStems
    gridStems{1}=temp;
end

%% Get mgrid info
[~, elecLabels, elecRgb]=mgrid2matlab(fsub,hem);

elecnames=cell(1,length(elecLabels));
for a=1:length(elecLabels),
    elecnames{a}=rmChar(elecLabels{a},'-');
end


%% Collect electrode neighbors for plotElecPial
pairs=cell(1,1);
uniStems=[];
uniStemsRgb=[];
ct=0;
for a=1:(length(elecLabels)-1)
    %     if isempty(findstr(lower(elecLabels{a}),'depth'))
    %         %not a depth electrode, find a neighbor
    id=find(elecLabels{a}=='-');
    elecnumA=str2num(elecLabels{a}(id+1:end));
    elecstemA=elecLabels{a}(1:id-1);
    
    if ~sum(ismember(lower(elecstemA),lower(gridStems)))
    %if ~sum(ismemberi(elecstemA,gridStems))
        % Ignore the grids
        if isempty(uniStems) || ~ismember(elecstemA,uniStems)
            uniStems{length(uniStems)+1}=elecstemA;
            uniStemsRgb(length(uniStems),1:3)=elecRgb(a,:);
        end
        
        id=find(elecLabels{a+1}=='-');
        elecnumB=str2num(elecLabels{a+1}(id+1:end));
        elecstemB=elecLabels{a+1}(1:id-1);
        if isempty(uniStems) || ~ismember(elecstemB,uniStems)
            uniStems{length(uniStems)+1}=elecstemB;
            uniStemsRgb(length(uniStems),1:3)=elecRgb(a+1,:);
        end
        
        if strcmp(elecstemA,elecstemB)
            ct=ct+1;
            pairs{ct,1}=elecnames{a};
            pairs{ct,2}=elecnames{a+1};
            pairs{ct,3}=elecRgb(a,:);
        end
    end
end

% Add grids
if ~isempty(gridStems)
    nGrid=length(gridStems);
    for a=1:nGrid,
        uniStems{length(uniStems)+1}=gridStems{a};
        % find the grid's color
        for b=1:length(elecLabels)
            if findstr(lower(gridStems{a}),lower(elecLabels{b}))
                uniStemsRgb(length(uniStems),1:3)=elecRgb(b,:);
                break;
            end
        end
        
        gridLines=derive_grid_lines(gridStems{a},gridDims(a,:));
        for b=1:size(gridLines,1)
            tempVec=gridLines{b,2};
            for c=2:length(tempVec)
                ct=ct+1;
                pairs{ct,1}=[gridLines{a,1} num2str(tempVec(c-1))];
                pairs{ct,2}=[gridLines{a,1} num2str(tempVec(c))];
                %pairs{ct,3}=[1 0 0]; % FIX!!
                pairs{ct,3}=uniStemsRgb(length(uniStems),:);
            end
        end
    end
end

% Format stems
nUni=length(uniStems);
for a=1:nUni,
     uniStems{a}=rmSubstring(uniStems{a},'depth');
end

%%
%for fLoop=1:1,
for fLoop=1:2,
    cfg=[];
    cfg.view=[hem 'omni'];
    cfg.figid=fLoop;
    cfg.eleccolors=elecRgb;
    cfg.elecnames=elecnames;
    cfg.plotcbar='n';
    cfg.ignore_depth_elec='n';
    cfg.pairs=pairs;
    %cfg.showlabels='y';
    if fLoop==2
        cfg.overlay_parcellation='DK';
    end
    %cfg.rotate3d='n';
    cfg.title=[];
    cfg_out=plotElecPial(fsub,cfg);
    
    % Add electrode legend
    hAx=axes('position',[.9 .01 .07 .98]);
    v=axis;
    dlt=(v(4)-v(3))/(nUni+1);
    for a=1:nUni,
        ht=text(v(2),v(3)+dlt*a,uniStems{a});
        set(ht,'fontsize',18,'color',uniStemsRgb(a,:), ...
            'horizontalalignment','right','fontweight','bold', ...
            'backgroundcolor','k');
    end
    set(hAx,'box','off','visible','off');
    
    set(gcf,'paperpositionmode','auto');
    if universalYes(printEm)
        if fLoop==1,
            figFname=sprintf('%s%sMgridElec',fsub,hem);
        else
            figFname=sprintf('%s%sMgridElecDK',fsub,hem);
        end
        fsubDir=getenv('SUBJECTS_DIR');
        print(fLoop,[fsubDir '/' fsub '/elec_recon/' figFname],'-djpeg');
    end
end

