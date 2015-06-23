function yangElecPjct(sub,hem)
%function yangElecPjct(sub,hem)
% 
% Corrects intracranial electrode locations for brain shift using the
% following method:
%
% Inputs:
%  sub - Freesurfer subject name (e.g., 'TWH001')
%  hem - ['lh' or 'rh'] Hemisphere in which electrodes were implanted.
%
% Outputs:
%  The following files are created in the elec_recon subfolder of the
%  Freesufer subject folder:
%    *_hem.PIAL - RAS coordinates snapped to pial surface
%    *_hem.PIALVOX - Voxel coordinates snapped to pial surface
%    *_hem.DURAL - RAS coordinates snapped to dural (i.e., smoothed pial)
%                  surface)
%    *_hem.electrodeNames - electrode names
%    localization_process_date.log - Record of command line output produced
%                                    when this function is run
% 
% In the above, *=Freesurfer subject name, hem='left' or 'right', and
% date=the date on which those files were generated
%
% Note, depth electrode coordinates are not affected by this function. They
% are kept the same as in the postimplant CT or MRI scan.
%
% Also note that you should double check grid labels as there may be a
% mismatch between NYU interpolation and mgrid file
%
% Author:
% David M. Groppe based on code from Hugh Wang
% Honeylab, University of Toronto
% June 2015

% Future work:
% -get mapping to MNI brain to work

% get the subject info
global global_fs_dir;
if ~isempty(global_fs_dir)
    fsDir=global_fs_dir;
else
    if ispc,
        error('Hey mon, if you be using Windows you need to be specifying global variable global_fs_dir.');
    else
        fsDir=getenv('SUBJECTS_DIR');
    end
end

subPath = sprintf('%s/%s',fsDir,sub);
elecReconPath=[subPath '/elec_recon/'];

if ~isdir(subPath)
    error('Freesurfer folder %s not found',subPath);
end

ii = strfind(subPath,'/');
fsDir = subPath(1:ii(end));


%% read the inital text file with postimplant electrode coordinates
% Note coordinates are in voxels (NOT RAS)
if strcmpi(hem,'lh') || strcmpi(hem,'l')
    hem='lh';
    postimpLocFname=sprintf('%s%sPostimpLocLeft.txt',elecReconPath,sub);
elseif strcmpi(hem,'rh') || strcmpi(hem,'r')
    hem='rh';
    postimpLocFname=sprintf('%s%sPostimpLocRight.txt',elecReconPath,sub);
else
   error('Illegal value of hem argument.'); 
end
if ~exist(postimpLocFname,'file')
    error('File %s does NOT exist.',postimpLocFname);
end


%%
% The code below migh be necessary for doing MNI mapping ?? DG commented
% them out
% [dep_img_file,dep_img_path] = uigetfile({'*.nii.gz';'*.nii';'*.img'},'Select the T1 pre-operation image: ',PathName);
% dep_img_path='/Applications/freesurfer/subjects/TWH008/elec_recon/';
% dep_img_file='T1.nii.gz';


if 0 % Not clear if this necessary
    % move other files to backup folder
    backup_dir = [PathName,'backup_',datestr(now,29)];
    if ~exist(backup_dir,'dir'), mkdir(backup_dir); end;
    all_files = dir(PathName);
    all_files_isdir = [all_files(:).isdir];
    all_files_name = {all_files(:).name};
    all_files_name = all_files_name(all_files_isdir==0);
    idx = ~cellfun(@isempty,strfind(all_files_name,'.log'))+...
        ~cellfun(@isempty,strfind(all_files_name,'pial_surf.mat'))+...
        ~cellfun(@isempty,strfind(all_files_name,'.ppt'))+...
        ~cellfun(@isempty,strfind(all_files_name,'.pptx'))+...
        ~cellfun(@isempty,strfind(all_files_name,FileName))+...
        ~cellfun(@isempty,strfind(all_files_name,dep_img_file))+...
        ~cellfun(@isempty,strfind(all_files_name,'missing')); % NYU settings: keep the NYxxx_missing_coor.txt
    all_files_name = all_files_name(~logical(idx));
    for i=1:length(all_files_name), movefile([PathName, all_files_name{i}],backup_dir,'f'); end;
end

% start diary
diary_file = [elecReconPath 'localization_process_',datestr(now,29),'.log'];
fprintf('Recording command line output in file: \n%s\n',diary_file);
diary on;
diary(diary_file)

fprintf('\n================================================================\n');
fprintf('Starting localization process for %s at %s\n',sub,datestr(now,31));
fprintf('Freesurfer Recon dir: %s\n',elecReconPath);
fprintf('Initial location text file: %s\n',postimpLocFname);

% read in preop T1
t1Fname=[elecReconPath '/T1.nii.gz'];
hdr = ntools_elec_load_nifti(t1Fname,1);
if ~isequal(hdr.pixdim(2),hdr.pixdim(3),hdr.pixdim(4))
    warning('T1 voxel mm dimensions not equal. Will affect the accuracy of distance calculation');
end
scale = mean(hdr.pixdim(2:4));

% read in postimp elec locations text/xls file
[~,~,ext] = fileparts(postimpLocFname);

if strcmpi(ext,'.txt')
    %fid = fopen([PathName, FileName]);
    fid = fopen(postimpLocFname);
    elec_all = textscan(fid,'%s %f %f %f %s','CommentStyle','%');
    elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4}),elec_all{5}];
    fclose(fid);
elseif strcmpi(ext,'.xls') || strcmpi(ext,'.xlsx')
    % I haven't tested this to make sure that it works. DG
    [~,~,elec_cell] = xlsread(postimpLocFname);
end

% split elecs by type
g = strmatch('G',upper(elec_cell(:,5)));
d = strmatch('D',upper(elec_cell(:,5)));
s = strmatch('S',upper(elec_cell(:,5)));

nManGrid=length(g);
nStrip=length(s);
nDepth=length(d);
nManElec=nDepth+nManGrid+nStrip;

fprintf('%d Depth electrodes manually marked\n',nDepth);
fprintf('%d Grid electrodes manually marked\n',nManGrid);
fprintf('%d Strip electrodes manually marked\n',nStrip);

% Convert all coordinates to RAS
postimpVox=zeros(nManElec,3);
for a=1:nManElec,
   for b=1:3,
      postimpVox(a,b)=elec_cell{a,b+1};
   end
end
VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];
postimpRas=(VOX2RAS*[postimpVox'; ones(1, nManElec)])';
%fprintf('RAS coordinates:\n');
for a=1:nManElec,
   for b=1:3,
      elec_cell{a,b+1}=postimpRas(a,b); 
   end
   %fprintf('%s: %.2f %.2f %.2f\n',elec_cell{a,1},elec_cell{a,2}, ...
   %    elec_cell{a,3},elec_cell{a,4});
end
elec_depth=elec_cell(d,:);
elec_depthCT=elec_depth;
elec_gridCT=elec_cell(g,:);
elec_stripCT=elec_cell(s,:);
elec_cell([g;d],:) = [];
nGrid=size(elec_gridCT,1);

%% Calculate the electrodes locations

% outer-brain surface check and create
ntools_elec_outer_brain(subPath); % If smoothed pial surface has not been created, I believe this function fails. DG

% calculate grids
if nManGrid
    % Figure out how many unique types of grids there are
    gridNames = regexp(elec_gridCT(:,1),'[A-Za-z]*[^\d*]','match');
    for i=1:length(gridNames)
        ini_gridNames(i) = gridNames{i};
    end
    gridNames = unique(ini_gridNames); %# of unique grids
    nGridType=length(gridNames);
    elec_grid=[];
    grid_stats=[];
    
    for a=1:nGridType,
        % For each grid identify the dimensions
        nRow=input(sprintf('How many rows does %s have? (default-> 8): ',gridNames{a}));
        if isempty(nRow),
           nRow=8; 
        end
        nCol=input(sprintf('How many columns does %s have? (default-> 8): ',gridNames{a}));
        if isempty(nCol),
            nCol=8;
        end
        
        % For each grid identify the corners in clockwise or counter clockwise
        % order
        defaultCorners=[1 nCol nCol*nRow nCol*nRow-nCol+1];
        corners=input(sprintf('Enter the electrode #''s of %s''s corners starting at 1 and going counterclockwise with square brackets (default-> [%d %d %d %d])', ...
            gridNames{a},defaultCorners(1),defaultCorners(2),defaultCorners(3),defaultCorners(4)));
        if isempty(corners),
            %corners=[1 8 64 57];
            corners=defaultCorners;
        end
        % Select out the corner electrodes to pass to Hugh's function:
        cornerIds=zeros(4,1);
        for b=1:4,
            cornerIds(b)=findstrInCell([gridNames{a} int2str(corners(b))],elec_gridCT(:,1),1);
        end
        
        %radius = input('\nInput the inter-electrode distance (mm) (Press Enter for default distance 10mm) : ');
        %if isempty(radius)
        %    radius = 10;
        %end
        radius=10;
        [elec_gridTemp, grid_statsTemp] = ntoolsElecCalcGrid(elec_gridCT(cornerIds,:),subPath,scale,radius,nRow,nCol);
        
        % Collect all temp grid coordinates
        elec_grid=[elec_grid; elec_gridTemp];
    end
end
%nGrid=size(elec_grid,1); % Right now, all grid elecs must be labeled in CT
%scan. If we change this so that only corners are labelled. We have to
%recount the number of grid elecs here. DG 
nElec=nDepth+nGrid+nStrip;

% calculate depth elecs
if 1
    fprintf('Localizing depths using manually marked locations.\n');
    fprintf('Locations are NOT corrected for brain shift.\n');
else
    fprintf('Localizing depths using the most extreme electrodes and interpolating the rest.\n');
    elec_depth = ntools_elec_calc_depth(ini_depth);
end

% calculate strip elecs
elec_strip = ntools_elec_calc_strip(elec_cell,subPath,hem);


%% Snap electrodes to the pial surface
pialRAS=zeros(nGrid+nStrip,3);
ct=0;
for a=1:nGrid,
    ct=ct+1;
    for b=1:3,
        pialRAS(ct,b)=elec_grid{a,b+1};
    end
end
for a=1:nStrip,
    ct=ct+1;
    for b=1:3,
        pialRAS(ct,b)=elec_strip{a,b+1};
    end
end
pialRAS = snap2surf(pialRAS,[subPath '/surf'],hem(1),'pial');

% add depth coords, which are not snapped to surface
for a=1:nDepth,
    ct=ct+1;
    for b=1:3,
        pialRAS(ct,b)=elec_depth{a,b+1};
    end
end


%% Save the electrodes locations into a text file
if strcmpi(hem(1),'l')
    longHem='left';
else
    longHem='right';
end

duralRAS=zeros(nElec,3);
ctRAS=duralRAS;
elecNames=cell(nElec,1);
ct=0;
for a=1:nGrid,
    ct=ct+1;
    elecNames{ct}=elec_grid{a,1};
    for b=1:3,
        duralRAS(ct,b)=elec_grid{a,b+1};
        tempId=findstrInCell(elecNames{ct},elec_gridCT,1);
        ctRAS(ct,b)=elec_gridCT{tempId,b+1};
    end
end
for a=1:nStrip,
    ct=ct+1;
    for b=1:3,
        duralRAS(ct,b)=elec_strip{a,b+1};
        ctRAS(ct,b)=elec_stripCT{a,b+1};
    end
    elecNames{ct}=elec_strip{a,1};
end
for a=1:nDepth,
    ct=ct+1;
    for b=1:3,
        duralRAS(ct,b)=elec_depth{a,b+1};
        ctRAS(ct,b)=elec_depthCT{a,b+1};
    end
    elecNames{ct}=elec_depth{a,1};
end

% RAS COORDINATES
% Dural
fnameDuralRAS = [elecReconPath sub '_' longHem '.DURAL'];
fprintf('Saving dural RAS electrode locations to: %s\n',fnameDuralRAS);
fid=fopen(fnameDuralRAS,'w');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'R A S\n');
for a=1:nElec,
   fprintf(fid,'%f %f %f\n',duralRAS(a,1),duralRAS(a,2),duralRAS(a,3)); 
end
fclose(fid);

% Pial
fnamePialRAS = [elecReconPath sub '_' longHem '.PIAL'];
fprintf('Saving pial RAS electrode locations to: %s\n',fnamePialRAS);
fid=fopen(fnamePialRAS,'w');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'R A S\n');
for a=1:nElec,
   fprintf(fid,'%f %f %f\n',pialRAS(a,1),pialRAS(a,2),pialRAS(a,3)); 
end
fclose(fid);

% CT (i.e., uncorrected for brain shift)
fnameCtRAS = [elecReconPath sub '_' longHem '.CT'];
fprintf('Saving CT RAS electrode locations to: %s\n',fnameCtRAS);
fid=fopen(fnameCtRAS,'w');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'R A S\n');
for a=1:nElec,
   fprintf(fid,'%f %f %f\n',ctRAS(a,1),ctRAS(a,2),ctRAS(a,3)); 
end
fclose(fid);

% Electrode names
fnameLabels = [elecReconPath sub '_' longHem '.electrodeNames'];
fprintf('Saving electrode labels to: %s\n',fnameLabels);
fid=fopen(fnameLabels,'w');
fprintf(fid,'%s\n',datestr(now));
for a=1:nElec,
   fprintf(fid,'%s\n',elecNames{a}); 
end
fclose(fid);

% VOX COORDINATES
RAS2VOX=inv(VOX2RAS);
duralVOX=(RAS2VOX*[duralRAS'; ones(1, nElec)])';
fnameDuralVOX = [elecReconPath sub '_' longHem '.DURALVOX'];
fprintf('Saving dural VOX electrode locations to: %s\n',fnameDuralVOX);
fid=fopen(fnameDuralVOX,'w');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'X Y Z\n');
for a=1:nElec,
    fprintf(fid,'%f %f %f\n',duralVOX(a,1),duralVOX(a,2),duralVOX(a,3));
end
fclose(fid);

pialVOX=(RAS2VOX*[pialRAS'; ones(1, nElec)])';
fnamePialVOX = [elecReconPath sub '_' longHem '.PIALVOX'];
fprintf('Saving pial VOX electrode locations to: %s\n',fnamePialVOX);
fid=fopen(fnamePialVOX,'w');
fprintf(fid,'%s\n',datestr(now));
fprintf(fid,'X Y Z\n');
for a=1:nElec,
    fprintf(fid,'%f %f %f\n',pialVOX(a,1),pialVOX(a,2),pialVOX(a,3));
end
fclose(fid);


%% save all into binary nifti image
if 0
    %fname_bin = [PathName,Sname,'_elec_bin_T1_' datestr(now,29), '.nii.gz'];
    fname_bin = [elecReconPath sub '_elec_bin_T1_' datestr(now,29), '.nii.gz'];
    fprintf('Saving electrode locations as binary nii file: %s\n',fname_bin);
    elec_vox = ntools_elec_savebin([x y z],hdr,fname_bin);
end

% transform into MNI space  ?? make work later
if 0
    elec_mni = ntools_elec_dartel_warp(fname_bin,[dep_img_path,dep_img_file]);
    %fname_mni = [PathName Sname '_coor_MNI_' datestr(now,29) '.txt'];
    fname_mni = [elecReconPath sub '_coor_MNI_' datestr(now,29) '.txt'];
    ntools_elec_savetxt(fname_mni,[name num2cell(elec_mni) label]);
end


%% Plot results to double check
shiftDist=sqrt( sum( (duralRAS-ctRAS).^2,2)); %units are mm

%rgb=zeros(nElec,3);
rgb=vals2Colormap(shiftDist,'justpos');
h_fig=figure;
%set(h_fig,'position',[104 285 1114 410]);
subplot(121);
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
subplot(122);
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

outFigFname=sprintf('%s/elec_recon/%s_%sShiftDist.jpg',subPath,sub,longHem);
print(gcf,'-djpeg',outFigFname);
outFigFname=sprintf('%s/elec_recon/%s_%sShiftDist',subPath,sub,longHem);
savefig(outFigFname);

%%
cfg=[];
cfg.view=[hem(1) 'omni'];
cfg.rotate3d='n';
cfg.eleccolors=shiftDist;
cfg.colorscale='justpos';
cfg.units='mm';
cfg.elecnames=elecNames;
%cfg.showlabels='y';
cfg.title=sprintf('%s: CT to Dural distance',sub);
cfg_out=plotElecPial(sub,cfg);

outFigFname=sprintf('%s/elec_recon/%s_%sShiftDistOnBrain.jpg',subPath,sub,longHem);
print(gcf,'-djpeg',outFigFname);

% close diary
fprintf('\nElectrodes Localization finished for %s',sub);
fprintf('\n================================================================\n');
diary off


end

