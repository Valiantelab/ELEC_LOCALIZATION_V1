function plotMgridOnPial(fsub,hem,printEm)
%function plotMgridOnPial(fsub,hem,printEm)
% 
% Makes romni or lomni plots on a gray pial and DK atlas pial of a
% patient's brain with electrodes colored according to mgrid file colors.
%
% Inputs:
%  fsub  - FreeSurfer subject directory
%  hem   - ['r' or 'l'] Hemisphere to plot
%  printEm - If nonzero jpgs of the figures will be made in the subject's
%            elec_recon folder
%
% Examples:
% plotMgridOnPial('TWH10','l',0);
%
%
% Author:
% David M. Groppe
% March, 2015
%

% Future Work:
%  -Add an option to make depths invisible

% History:
% May 2015-Now electrode pairs are derived directly from mgrid file and
% disabled electrodes are not shown.

if nargin<4
    gridStems=[];
elseif ~iscell(gridStems)
    %convert string to a singleton cell array
    temp=gridStems;
    clear gridStems
    gridStems{1}=temp;
end

%% Get mgrid info
%[~, elecLabels, elecRgb]=mgrid2matlab(fsub,hem);
[~, elecLabels, elecRgb, elecPairs, elecPresent]=mgrid2matlab(fsub,hem);

elecnames=cell(1,length(elecLabels));
for a=1:length(elecLabels),
    elecnames{a}=rmChar(elecLabels{a},'-');
end

pairPresent=zeros(size(elecPairs,1),1);
for a=1:size(elecPairs,1),
   elecPairs{a,1}=rmChar(elecPairs{a,1},'-');
   elecPairs{a,2}=rmChar(elecPairs{a,2},'-');
   elecId1=findstr_in_cell(elecPairs{a,1},elecnames,1);
   elecId2=findstr_in_cell(elecPairs{a,2},elecnames,1);
   pairPresent(a)=elecPresent(elecId1)*elecPresent(elecId2);
end

%% Collect electrode neighbors for plotElecPial
uniStems=[];
uniStemsRgb=[];
ct=0;
for a=1:length(elecLabels),
    %     if isempty(findstr(lower(elecLabels{a}),'depth'))
    %         %not a depth electrode, find a neighbor
    id=find(elecLabels{a}=='-');
    elecstem=elecLabels{a}(1:id-1);
    
    if ~ismember(elecstem,uniStems)
        ct=ct+1;
        uniStems{ct}=elecstem;
        uniStemsRgb(ct,1:3)=elecRgb(a,:);
    end

end


% Format stems
nUni=length(uniStems);
for a=1:nUni,
     uniStems{a}=rmSubstring(uniStems{a},'depth');
end


%%
%for fLoop=1:1, % ??
for fLoop=1:2,
    cfg=[];
    cfg.view=[hem 'omni'];
    cfg.figid=fLoop;
    cfg.eleccolors=elecRgb(find(elecPresent),:);
    cfg.elecnames=elecnames(find(elecPresent));
    cfg.plotcbar='n';
    cfg.ignore_depth_elec='n';
    cfg.pairs=elecPairs(find(pairPresent),:);
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

