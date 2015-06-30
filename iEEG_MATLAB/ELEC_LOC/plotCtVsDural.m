function [figH, axH]=plotCtVsDural(duralRAS,ctRAS,elecNames,sub,hem,subPath,printEm,plotPial)
%function [figH, axH]=plotCtVsDural(duralRAS,ctRAS,elecNames,sub,hem,subPath,printEm,plotPial)
%
% This function creates two plots to illustrate the effect of brain shift
% correction:
%  1: A plot of the distance between each electrode's pre and post
%  correction positions and a 3D scatter plot of these positions
%
%  2: Post-correction positions overlayed on pial surface and color coded
%  to represent distance between pre and post correction locations
%
% This function is called pial interpStripElec.m and yangElecPjct.m
%
% Author: David Groppe
% Honeylab, Univ. of Toronto
% June 2015

%% Plot results to double check
shiftDist=sqrt( sum( (duralRAS-ctRAS).^2,2)); %units are mm

%rgb=zeros(nElec,3);
rgb=vals2Colormap(shiftDist,'justpos');
figH(1)=figure;
%set(figH(1),'position',[104 285 1114 410]);
axH(1)=subplot(121);
plot(shiftDist,'.-'); hold on;
last_type=[];
marker='o';
non_depth=ones(1,length(shiftDist));
for a=1:length(shiftDist),
    if ~isempty(findstr('depth',elecNames{a}))
        non_depth(a)=0;
        marker='s';
    else
        marker='o';
    end
    %    ids=find(elecNames{a}=='_');
    %    this_type=elecNames{a}(1:(ids(1)-1));
    %    if ~strcmpi(this_type,last_type)
    %        last_type=this_type;
    %        if marker=='o';
    %            marker='s';
    %        else
    %            marker='o';
    %        end
    %    end
    h=plot(a,shiftDist(a),marker);
    set(h,'color',rgb(a,:));
    clickText(h,[num2str(a) ': ' rmChar(elecNames{a},'_')]);
    %clickText(h,[num2str(a) ': ' rm_substring(elecNames{a},'_')]);
    xlabel('Channel');
    ylabel('Snapped vs. Nonsnapped Distance (mm)');
end
non_depth_ids=find(non_depth);
title(sprintf('%s Median=%.1f, SIQR=%.1f (depths ignored)',sub,median(shiftDist(non_depth_ids)), ...
    iqr(shiftDist(non_depth_ids))/2));
v=axis;
axis([1 length(shiftDist) v(3:4)]);

%3D plot of and pre vs. post shift correction locations
axH(2)=subplot(122);
for a=1:length(shiftDist),
    h=plot3(duralRAS(a,1),duralRAS(a,2),duralRAS(a,3),'r.'); hold on;
    clickText(h,rmChar(elecNames{a},'_'));
    h=plot3(ctRAS(a,1),ctRAS(a,2),ctRAS(a,3),'bo');
    %clickText(h,rm_substring(labels{a},'_'));
    clickText(h,rmChar(elecNames{a},'_'));
    plot3([duralRAS(a,1) ctRAS(a,1)],[duralRAS(a,2) ctRAS(a,2)],[duralRAS(a,3) ctRAS(a,3)],'k-');
end
if hem(1)=='r',
    %view(90,0);
    view([100 16]);
else
    %view(270,0);
    view([-76 12]);
end
axis tight;
axis square;
title('Red=Postcorrection, Blue=Precorrection');
xlabel('Left- Right+');
ylabel('Pos- Ant+');
zlabel('Inf- Sup+');

% Plot shift distances on pial surface
if universalYes(plotPial)
    figH(2)=figure;
    cfg=[];
    cfg.view=[hem(1) 'omni'];
    cfg.rotate3d='n';
    cfg.figid=figH(2);
    cfg.eleccolors=shiftDist;
    cfg.colorscale='justpos';
    cfg.units='mm';
    cfg.elecnames=elecNames;
    %cfg.showlabels='y';
    cfg.title=sprintf('%s: CT to Dural distance',sub);
    cfg_out=plotElecPial(sub,cfg);
    
    if universalYes(printEm)
        if strcmpi(hem,'r') || strcmpi(hem,'rh')
            longHem='right';
        elseif strcmpi(hem,'l') || strcmpi(hem,'lh')
            longHem='left';
        else
            error('Invalid value for "hem" argument.');
        end
        outFigFname=sprintf('%s/elec_recon/%s_%sShiftDist.jpg',subPath,sub,longHem);
        print(gcf,'-djpeg',outFigFname);
        outFigFname=sprintf('%s/elec_recon/%s_%sShiftDist',subPath,sub,longHem);
        savefig(outFigFname);
        outFigFname=sprintf('%s/elec_recon/%s_%sShiftDistOnBrain.jpg',subPath,sub,longHem);
        print(gcf,'-djpeg',outFigFname);
    end
end
