function [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(mgridFname,hem)
%function [elecMatrix, elecLabels, elecRgb]=mgrid2matlab(mgridFname,hem)
%
% Optional Inputs:
%  mgridFname - The filename AND path of the mgrid file you wish to import or the
%               freesurfer name of the subject. If not specified, an 
%               spm_select gui pops up for you to save the file
%  hem        - ['l' or 'r'] The hemisphere of the mgrid file you wish to
%               load. Necessary (and only used) if mgridFname is a freesurfer name.
%
% Output:
%  elecMatrix - 3D matrix of electrode coordinates. Column 1 is R->L (i.e.,
%               1 is the right-most plane of voxels). Column 2 is in S->I
%               (i.e., 1 is the most superior plane of voxels). Column 3 is
%               A->P (i.e., 1 is the most anterior plane of voxels).
%  elecLabels - Cell array of electrode names corresponding to each row of
%               elecMatrix
%  elecRgb    - Matrix of RGB colors for each electrode. The nth row of
%               elecRgb corresponds to the nth element of elecLabels.
%
%  Examples:
%   Load with a GUI:
%   >>[elecMatrix, elecLabels, elecRgb]=mgrid2matlab();
%
%   Load using FreeSurfer storage conventions:
%   >>[elecMatrix, elecLabels, elecRgb]=mgrid2matlab('TWH11','l');
%
%   Load using an full mgrid filename:
%   >>[elecMatrix, elecLabels, elecRgb]=mgrid2matlab('/Applications/freesurfer/subjects/TWH11/elec_recon/TWH11_left.mgrid');
%
%
% Authors: Saba Shahab & David M. Groppe
% Feb. 2015
% Feinstein Institute for Medical Research/Univ. of Toronto

% Future work:
% -get this function to return electrode neighbors too?
% -check to see if it works for invisible electrodes?

if nargin<1
    [dataMGRID.elecName]=spm_select(1,'mgrid','select the .mgrid file generated by BioImageSuite');
else
    fsdir=getenv('SUBJECTS_DIR');
    if strcmpi(hem,'L')
        files=dir(sprintf('%s/%s/elec_recon/*left.mgrid',fsdir,mgridFname));
    elseif strcmpi(hem,'R')
        files=dir(sprintf('%s/%s/elec_recon/*right.mgrid',fsdir,mgridFname));
    else
        error('Invalid value for "hem" parameter.');
    end
    if length(files)>1
        error('Multiple mgrid files exist for this hemisphere.');
    elseif isempty(files)
        error('No mgrid files exist for this hemisphere.');
    end
    mgridFname=sprintf(sprintf('%s/%s/elec_recon/%s',fsdir,mgridFname,files(1).name));
end
fid = fopen(mgridFname, 'r');

elecMatrix = [];
elecLabels=[];
elecRgb=[];
ct=0;
crntLabel=[];
crntCt=0;
crntColor=[1 1 1];
while feof(fid) == 0
    line = fgetl(fid);
    if strcmpi(line, '#Position')
        line = fgetl(fid);
        ct=ct+1;
        elecMatrix=[elecMatrix; str2num(line)];
        elecLabels{ct}=[crntLabel '-' num2str(crntCt)];
        elecRgb(ct,1:3)=crntColor;
        crntCt=crntCt-1;
    elseif strcmpi(line,'#Description')
        crntLabel = fgetl(fid);
    elseif strcmpi(line,'#Dimensions')
        dim = str2num(fgetl(fid));
        crntCt=dim(1)*dim(2);
    elseif strcmpi(line,'#Color')
        crntColor=str2num(fgetl(fid));
    end
    %elseif strncmpi(line,'# Electrode Grid')
end
fid = fclose(fid);
elecMatrix=elecMatrix+1; % BioImageSuite indexes the first voxel as [0 0 0]