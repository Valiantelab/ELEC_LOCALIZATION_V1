function [inf_coor, labels]=pvox2InfBrain(subj,hem,cfg)
%function [avg_coords, avg_vids, sub_vids, labels]=pvox2InfBrain(subj,hem,cfg)
%

% This function takes the "pial" coordinates (snapped to pial surface) and:
% 1. finds the closest vertex on the pial surface
% 2. maps that vertex to the inflated brain surface for that subject
%
% Inputs:
%   subj = FreeSurfer subject name
%   hem  ='l' or 'r'
%
% Optional Inputs: passed as fields in a configuration structure
%   eleccoord = N-by-3 numeric array with electrode RAS coordinates. {default:
%               not used; the function looks into the subject's Freesurfer
%               folder for electrode coordinate file instead}
%   elecnames = cell array of strings with electrode names, corresponding
%               to the rows of eleccoord. {default: not used; the function
%               looks into the subject's Freesurfer folder for electrode
%               name file instead}
%   fsurfsubdir = path to the Freesurfer subject directory. Necessary if
%                 running MATLAB on Windows. {default: taken from shell}
%
% Outputs:
%   inf_coords = Electrode coordinates on FreeSurfer inflated pial surface
%   labels     = Channel names (from *.PIAL file)
%
%
% Author: 
% David Groppe
% Mehtalab
% April, 2013
%

% History:
% 2015-6 Made compatible with new Yang method for brain shift correction: DG

% parse input parameters in cfg structure and set defaults
if  ~isfield(cfg,'eleccoord'),      eleccoord = []; else    eleccoord = cfg.eleccoord;      end
if  ~isfield(cfg,'elecnames'),      elecnames = []; else    elecnames = cfg.elecnames;      end
if  ~isfield(cfg,'fsurfsubdir'),    fs_dir = [];    else    fs_dir = cfg.fsurfsubdir;       end

hem=lower(hem);

% Get location of FreeSurfer directories
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
sub_dir=[fs_dir '/' subj];

%% Read Sub Pial Surf
fname=[sub_dir '/surf/' hem 'h.pial'];
pial=readSurfHelper(fname);

%% get electrode coordinates
if isempty(eleccoord) % no electrode coordinates have been passed in the function call: use the original code looking for .PIAL files
    f=dir([sub_dir '/elec_recon/*.PIAL']);
    files4bothsides=0;
    if length(f)>1,
        if hem(1)=='l',
            f=dir([sub_dir '/elec_recon/*left.PIAL']);
            files4bothsides=1;
        else
            f=dir([sub_dir '/elec_recon/*right.PIAL']);
            files4bothsides=1;
        end
    end
    if length(f)>1,
        error('Too many possible PIAL files.  I do not know which to use.');
    end
    elec_coord=[sub_dir '/elec_recon/' f(1).name];
    id=find(f(1).name=='.');
    label_fname=[sub_dir '/elec_recon/' f(1).name(1:id) 'electrodeNames'];
    labels=textread(label_fname,'%s');
    labels=labels(3:end); %get rid of data header

    tempCsv=csv2cell(elec_coord,' ',2);
    RAS_coor=zeros(size(tempCsv,1),1);
    for csvLoopA=1:size(tempCsv,1),
        for csvLoopB=1:3,
            RAS_coor(csvLoopA,csvLoopB)=str2num(tempCsv{csvLoopA,csvLoopB});
        end
    end
    if length(labels)~=size(RAS_coor,1)
        error('# of electrode labels does NOT equal # of coordinates in PIAL file.');
    end
else % numeric electrode coordinates have been passed in the function call
    RAS_coor=eleccoord;
    labels=elecnames;
end


%% Get vertices for electrodes on subject's pial surface
n_chan=size(RAS_coor,1);
n_pial_vert=size(pial.vert,1);
sub_vids=zeros(1,n_chan);
for a=1:n_chan,
    df=pial.vert-repmat(RAS_coor(a,:),n_pial_vert,1);
    dst=sum(abs(df),2);
    [dummy sub_vids(a)]=min(dst);
end

%% Read Sub Inflated Surf
fname=[sub_dir '/surf/' hem 'h.inflated'];
inflated=readSurfHelper(fname);

%% Get electrode coordinates on inflated brain
inf_coor=inflated.vert(sub_vids,:);





