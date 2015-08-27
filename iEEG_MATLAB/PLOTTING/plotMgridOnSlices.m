function plotMgridOnSlices(fsSub,mgridFname,cfg)
% function plotMgridOnSlices(fsSub,mgridFname,cfg)
%
% Creates a figure illustrating the location of each electrode in an mgrid
% file in a sagittal, coronal, and axial slice and indicates which part of
% the brain it is in.
%
% Required Inputs:
%  fsSub - Patient's freesurfer directory name
%  mgridFname - mgrid filename and path
%
% Optional cfg parameters:
%  fullTitle - If 1, the mgrid and mri voxel coordinates are displayed in
%              the figure title along with the electrode name and anatomical 
%              location. {default: 0}
%  cntrst    - 0< number <=1 The lower this number the lower the brighter
%              the image (i.e., the lower the voxel value corresponding to 
%              white). {default: 0.5}
%  pauseOn   - If 1, Matlab pauses after each figure is made and waits for
%              a keypress. {default: 0}
%  printFigs - If 1, each figure is output to an eps file. {default: 0}
%  depthsOnly - If 1, figures are only made for depth electrodes. Depth 
%               electrodes must have the string 'depth' in their name. 
%               {default: 1}
%
% Examples:
%  %Specify mgrid file and do NOT print
%  mgridFname='/Applications/freesurfer/subjects/NiAs/elec_recon/NiAs_bi.mgrid';
%  cfg=[];
%  plotMgridOnSlices('NiAs',mgridFname,cfg);
%
%  %Use FreeSurfer file structure and print
%  cfg=[];
%  cfg.printFigs=1;
%  plotMgridOnSlices('TWH11','l',cfg);
%
%
% Author: David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto

% Future work:
% Add option for fsurf anatomy colors?

if ~isfield(cfg,'fullTitle'),    fullTitle=0;          else fullTitle=cfg.fullTitle; end
if ~isfield(cfg,'cntrst'),    cntrst=.5;          else cntrst=cfg.cntrst; end
if ~isfield(cfg,'pauseOn'),    pauseOn=0;          else pauseOn=cfg.pauseOn; end
if ~isfield(cfg,'printFigs'),    printFigs=0;          else printFigs=cfg.printFigs; end
if ~isfield(cfg,'depthsOnly'),    depthsOnly=1;          else depthsOnly=cfg.depthsOnly; end
checkCfg(cfg,'plotMgridOnSlices.m');

% Load MRI
global global_fs_dir;
if ~isempty(global_fs_dir)
    fsdir=global_fs_dir;
else
    if ispc,
        error('Hey mon, if you be using Windows you need to be specifying global variable global_fs_dir.');
    else
        fsdir=getenv('SUBJECTS_DIR');
    end
end

mriFname=[fsdir '/' fsSub '/mri/brainmask.mgz'];
if ~exist(mriFname,'file')
   error('File %s not found.',mriFname); 
end
mri=MRIread(mriFname);
%mri.vol is ILA (i.e., S->I, R->L, P->A)
mx=max(max(max(mri.vol)))*cntrst;
mn=min(min(min(mri.vol)));
sVol=size(mri.vol);

% Load mgrid
if strcmpi(mgridFname,'l') || strcmpi(mgridFname,'r')
    [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(fsSub,mgridFname);
else
    [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(mgridFname); % mgrid coords are LIP
end
nElec=length(elecLabels);
elecMatrix=round(elecMatrix);
xyz=zeros(size(elecMatrix));
xyz(:,1)=elecMatrix(:,2);
xyz(:,2)=elecMatrix(:,1);
xyz(:,3)=sVol(3)-elecMatrix(:,3);


for elecId=1:nElec,
    if universalYes(depthsOnly)
        plotEm=~isempty(findstr('depth',lower(elecLabels{elecId})));
    else
        plotEm=1;
    end
    
    if plotEm
        figId=figure();
        set(figId,'position',[78 551 960 346],'paperpositionmode','auto');
        
        hm=zeros(1,3);
        figure(figId); clf;
        colormap gray;
        %subplot(131);
        wdth=.35;
        wDelt=.33;
        xStart=-.005;
        yStart=.03;
        ht=.9;
        axes('position',[xStart yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,xyz(elecId,2),:)),[mn mx]);
        axis square;
        set(gca,'xdir','reverse');
        hold on;
        hm(1)=plot(xyz(elecId,3),xyz(elecId,1),'r.');
        set(hm(1),'color',elecRgb(elecId,:));
        %find image limits
        mxX=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],2);
        mxY=max(squeeze(mri.vol(:,xyz(elecId,2),:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(1)/2),find(mxY==0)));
        limYb=min(intersect((sVol(1)/2:sVol(1)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        axis([tempMin tempMax tempMin tempMax]);
        set(gca,'xtick',[],'ytick',[]);
        
        %subplot(132);
        axes('position',[xStart+wDelt yStart wdth ht]);
        imagesc(squeeze(mri.vol(xyz(elecId,1),:,:)),[mn mx]);
        axis square;
        hold on;
        hm(2)=plot(xyz(elecId,3),xyz(elecId,2),'r.');
        set(hm(2),'color',elecRgb(elecId,:));
        %find image limits
        mxX=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],2);
        mxY=max(squeeze(mri.vol(xyz(elecId,1),:,:)),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        axis([tempMin tempMax tempMin tempMax]);
        set(gca,'xtick',[],'ytick',[],'xdir','reverse');
        
        
        %subplot(133);
        axes('position',[xStart+wDelt*2 yStart wdth ht]);
        imagesc(squeeze(mri.vol(:,:,xyz(elecId,3))),[mn mx]);
        axis square;
        hold on;
        hm(3)=plot(xyz(elecId,2),xyz(elecId,1),'r.');
        set(hm(3),'color',elecRgb(elecId,:));
        %find image limits
        mxX=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],2);
        mxY=max(squeeze(mri.vol(:,:,xyz(elecId,3))),[],1);
        limXa=max(intersect(1:(sVol(3)/2),find(mxX==0)));
        limXb=min(intersect((sVol(3)/2:sVol(3)),find(mxX==0)));
        limYa=max(intersect(1:(sVol(2)/2),find(mxY==0)));
        limYb=min(intersect((sVol(2)/2:sVol(2)),find(mxY==0)));
        %keep image square
        tempMin=min([limXa limYa]);
        tempMax=max([limXb limYb]);
        axis([tempMin tempMax tempMin tempMax]);
        set(gca,'xtick',[],'ytick',[]);
        
        anatLabel=vox2seg(xyz(elecId,:),fsSub);
  
        % Remove string "depth" from electrode label
        formattedLabel=rmSubstring(elecLabels{elecId},'depth',0);
        formattedLabel=rmChar(formattedLabel,'_');
        
        if universalYes(fullTitle)
            ht=textsc2014([formattedLabel '; mgrid coords(' num2str(elecMatrix(elecId,:)-1) '); fsurf coords(' num2str(xyz(elecId,:)) '); ' anatLabel], ...
                'title');
            set(ht,'fontsize',14,'fontweight','bold');
        else
            ht=textsc2014([formattedLabel '; Anatomical Location: ' anatLabel], ...
                'title');
            set(ht,'fontsize',16,'fontweight','bold');
        end
        set(ht,'position',[.5 .97 0]);
        
        if universalYes(printFigs)
            for a=1:3,
                set(hm(a),'markersize',14);
            end
            %print(fLoop,[fsdir '/' fsub '/elec_recon/' figFname],'-djpeg');
            %figFname=sprintf('%s_%sSlices',fsSub,elecLabels{elecId});
            figFname=sprintf('%s/%s/elec_recon/%s_%sSlices',fsdir,fsSub,fsSub,elecLabels{elecId});
            fprintf('Exporting figure to %s\n',figFname);
            %print(figId,figFname,'-depsc');
            print(figId,figFname,'-djpeg');
        end
        
        if universalYes(pauseOn)
            fprintf('Paused. Press any key for next electrode.\n');
            pause;
        end
    end
end
fprintf('Done showing all electrodes.\n');