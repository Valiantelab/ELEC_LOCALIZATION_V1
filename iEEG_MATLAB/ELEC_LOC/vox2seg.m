function anatLabel=vox2aseg(coordILA,fsSub)
%function anatLabel=vox2aseg(coordILA,fsSub)
%
% Inputs:
%  coordILA - 3D vector indicating the coordinates of a voxel in a
%             FreeSurfer MRI. First coordinate is Sup->Inf (i.e., 1=most
%             superior slice). Second coordinate is R->L (i.e., 1=the
%             rightmost slice). Third coordinate is P->A (i.e., 1=the
%             most posterior slice).
%  fsSub    - The FreeSurfer folder for the patient
%
% Output:
%  anatLabel - The part of the brain that voxel belongs to (e.g., Left
%              Hippocampus) according to aparc+aseg.mgz
%
%  Example:
%  >>coordILA=[147 144 115];
%  >>anatLabel=vox2aseg(coordILA,'NiAs')
%
% Author: David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto

% Load aseg volume
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
%asegFname=[fsdir '/' fsSub '/mri/aseg.mgz'];
asegFname=[fsdir '/' fsSub '/mri/aparc+aseg.mgz'];
%asegFname=[fsdir '/' fsSub '/mri/aparc.a2009s+aseg.mgz'];
if ~exist(asegFname,'file')
   error('File %s not found.',asegFname); 
end

aseg=MRIread(asegFname);

%% Load table
pathstr = fileparts(which('mgrid2matlab'));
inFile=[pathstr '/FreeSurferColorLUTnoFormat.txt'];
if ~exist(inFile,'file')
    error('Could not find file %s',inFile);
end
fid=fopen(inFile,'r');
%fid=fopen('/Applications/freesurfer/FreeSurferColorLUTnoFormat.txt','r');
tbl=textscan(fid,'%d%s%d%d%d%d');
fclose(fid);

%% Find anatomical region corresponding to voxel
id=find(aseg.vol(coordILA(1),coordILA(2),coordILA(3))==tbl{1});
id=min(id);
anatLabel=tbl{2}{id};
