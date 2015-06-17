function [avg_coords, avg_vids, sub_vids, labels]=pvox2AvgBrain(subj,hem,cfg)
%function [avg_coords, avg_vids, sub_vids, labels]=pvox2AvgBrain(subj,hem,cfg)
%
% This function takes RAS "pial" coordinates (snapped to pial surface)
% and maps it to the corresponding location on the pial
% surface of FreeSurfer's average brain.
%
% Inputs:
%   subj = FreeSurfer subject name
%   hem  ='l' or 'r'
%
% Optional Inputs: passed as fields in a configuration structure
%   plotem = 1 or 0.  If nonzero, a figure is created illustrating
%            electrode locations on subject and average pial surface.
%            Click on electrodes to see names. Depth electrodes are not
%            shown.  Code assumes the string 'depth' is part of the name of
%            all depth electrodes. {default: 1}
%   eleccoord = N-by-3 numeric array with RAS electrode coordinates. {default:
%               not used; the function looks into the subject's Freesurfer
%               folder for electrode coordinate file instead}
%   elecnames = cell array of strings with electrode names, corresponding
%               to the rows of eleccoord. {default: not used; the function
%               looks into the subject's Freesurfer folder for electrode
%               name file instead}
%   fsurfsubdir = path to the Freesurfer subject directory. Necessary if
%                 running MATLAB on Windows. {default: taken from shell}
%   rmdepths = 1 or 0. If nonzero, depth electrodes are ignored. Depth
%              elecrodes are indicated by having the string "depth" in 
%              their name. {default: 1}
%
% Outputs:
%   avg_coords = Electrode coordinates on FreeSurfer avg brain pial surface
%                (RAS coordinates)
%   avg_vids   = Subject pial surface vertices corresponding to each electrode
%   sub_vids   = Average pial surface vertices corresponding to each electrode
%   labels     = Channel names (from PIAL file)
%
%
% Author:
% David Groppe
% Mehtalab
% March, 2012
%

% History
% 2015-6 Made compatible with Yang brain shift correction algorithm

% parse input parameters in cfg structure and set defaults
if  ~isfield(cfg,'plotem'),         plotem = 1;     else    plotem = cfg.plotem;            end
if  ~isfield(cfg,'eleccoord'),      eleccoord = []; else    eleccoord = cfg.eleccoord;      end
if  ~isfield(cfg,'elecnames'),      elecnames = []; else    elecnames = cfg.elecnames;      end
if  ~isfield(cfg,'fsurfsubdir'),    fs_dir = [];    else    fs_dir = cfg.fsurfsubdir;       end
if  ~isfield(cfg,'rmdepths'),         rmdepths = 1;     else    rmdepths = cfg.rmdepths;            end

hem=lower(hem);

global global_fs_dir;
if isempty(fs_dir)
    if ~isempty(global_fs_dir)
        fs_dir=global_fs_dir;
    else
        if ispc,
            error('Hey mon, if you be using Windows you need to be specifying global_fs_dir.');
        else
            fs_dir=getenv('SUBJECTS_DIR');
        end
    end
end
avg_dir=[fs_dir '/' 'fsaverage'];
sub_dir=[fs_dir '/' subj];

fname=[sub_dir '/surf/' hem 'h.pial'];
pial=readSurfHelper(fname);

if isempty(eleccoord) % no electrode coordinates have been passed in the function call: use the original code looking for .PIAL files
    f=dir([sub_dir '/elec_recon/*.PIAL']);
    files4bothsides=0;
    if length(f)>1,
        if strcmp(hem,'r'),
            f=dir([sub_dir '/elec_recon/*right.PIAL']);
            files4bothsides=1;
        else
            f=dir([sub_dir '/elec_recon/*left.PIAL']);
            files4bothsides=1;
        end
    end
    if length(f)>1,
        error('Too many possible PIAL files. I do not know which to use.');
    end
    elec_coord_fname=[sub_dir '/elec_recon/' f(1).name];
    tempCsv=csv2Cell(elec_coord_fname,' ',2);
    RAS_coor=zeros(size(tempCsv,1),1);
    for csvLoopA=1:size(tempCsv,1),
        for csvLoopB=1:3,
            RAS_coor(csvLoopA,csvLoopB)=str2num(tempCsv{csvLoopA,csvLoopB});
        end
    end
    
    id=find(f(1).name=='.');
    label_fname=[sub_dir '/elec_recon/' f(1).name(1:id) 'electrodeNames'];
    labels=textread(label_fname,'%s');
    labels=labels(3:end); %get rid of data header
    
else % numeric electrode coordinates have been passed in the function call
    RAS_coor=eleccoord;
    labels=elecnames;
    if size(RAS_coor,1)~=length(labels)
        error('elecnames and eleccoord need to have the same number of electrodes.');
    end
end

n_chan=size(RAS_coor,1);
n_pial_vert=size(pial.vert,1);
sub_vids=zeros(1,n_chan);
for a=1:n_chan,
    df=pial.vert-repmat(RAS_coor(a,:),n_pial_vert,1);
    dst=sum(abs(df),2);
    [dummy sub_vids(a)]=min(dst);
end

fname=[sub_dir '/surf/' hem 'h.sphere.reg'];
sph=readSurfHelper(fname);

fname=[avg_dir '/surf/' hem 'h.sphere.reg'];
avg_sph=readSurfHelper(fname);
n_avg_vert=length(avg_sph.vert);
avg_vids=zeros(1,n_chan);
for b=1:n_chan,
    df=avg_sph.vert-repmat(sph.vert(sub_vids(b),:),n_avg_vert,1);
    dst=sum(abs(df),2);
    [dummy avg_vids(b)]=min(dst);
end


fname=[avg_dir '/surf/' hem 'h.pial'];
avg_pial=readSurfHelper(fname);
if min(min(avg_pial.tri))<1
    avg_pial.tri=avg_pial.tri+1; %sometimes this is needed sometimes not. no comprendo. DG ??
end
avg_coords=zeros(n_chan,3);
for a=1:n_chan,
    avg_coords(a,:)=avg_pial.vert(avg_vids(a),:);
end

% Remove depths
if universalYes(rmdepths)
    useChans=logical(false(n_chan,1));
    for a=1:n_chan
        if isempty(findstr('depth',lower(labels{a})))
            useChans(a)=1;
        end
    end
    avg_coords=avg_coords(useChans,:);
    avg_vids=avg_vids(useChans);
    sub_vids=sub_vids(useChans);
    labels=labels(useChans);
    n_chan=length(labels);
end

if universalYes(plotem)
    nice_labs=formatElecNames(labels);
    
    h_fig=figure;
    set(h_fig,'position',[360 335 829 360]);
    subplot(1,2,1);
    map=[1 1 1]*.7;
    tripatchDG(avg_pial,h_fig,map);
    shading interp; lighting gouraud; material dull; axis off, hold on
    if strcmp(hem,'r')
        l=light('Position',[1 0 0]);
        view(90,0);
    else
        l=light('Position',[-1 0 0]);
        view(270,0);
    end
    for a=1:n_chan
        h=plot3(avg_coords(a,1),avg_coords(a,2),avg_coords(a,3),'r.');
        clickText(h,nice_labs{a});
        set(h,'markersize',20);
    end
    rotate3d off;
    
    subplot(1,2,2);
    if min(min(pial.tri))<1
        pial.tri=pial.tri+1; %sometimes this is needed sometimes not. no comprendo. DG ??
    end
    tripatchDG(pial,h_fig,map);
    shading interp; lighting gouraud; material dull; axis off, hold on
    if strcmp(hem,'r')
        l=light('Position',[1 0 0]);
        view(90,0);
    else
        l=light('Position',[-1 0 0]);
        view(270,0);
    end
    for a=1:n_chan
        d=sub_vids(a);
        h=plot3(pial.vert(d,1),pial.vert(d,2),pial.vert(d,3),'r.');
        clickText(h,nice_labs{a});
        set(h,'markersize',20);
    end
    rotate3d off;
    set(gcf,'name',subj);
    
    h_fig=figure;
    set(h_fig,'position',[360 335 829 360]);
    subplot(1,2,1);
    map=[1 1 1]*.7;
    tripatchDG(avg_pial,h_fig,map);
    shading interp; lighting gouraud; material dull; axis off, hold on
    if strcmp(hem,'r')
        l=light('Position',[-1 0 0]);
        view(270,0)
    else
        l=light('Position',[1 0 0]);
        view(90,0);
    end
    for a=1:n_chan
        h=plot3(avg_coords(a,1),avg_coords(a,2),avg_coords(a,3),'r.');
        clickText(h,[nice_labs{a} sprintf(' %.3f %.3f %.3f',avg_coords(a,1), ...
            avg_coords(a,2),avg_coords(a,3))]);
        set(h,'markersize',20);
    end
    rotate3d off;
    
    subplot(1,2,2);
    if min(min(pial.tri))<1
        pial.tri=pial.tri+1; %sometimes this is needed sometimes not. no comprendo. DG ??
    end
    tripatchDG(pial,h_fig,map);
    shading interp; lighting gouraud; material dull; axis off, hold on
    if strcmp(hem,'r')
        l=light('Position',[-1 0 0]);
        view(270,0)
    else
        l=light('Position',[1 0 0]);
        view(90,0);
    end
    for a=1:n_chan
        d=sub_vids(a);
        h=plot3(pial.vert(d,1),pial.vert(d,2),pial.vert(d,3),'r.');
        clickText(h,nice_labs{a});
        set(h,'markersize',20);
    end
    rotate3d off;
    set(gcf,'name',subj);
    
    drawnow;
end
%avg_coords=[VOX2RAS\[avg_coords ones(n_chan,1)]']';
%avg_coords=avg_coords(:,1:3);