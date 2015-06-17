function makeIniLocTxtFile(sub,hem)
%function makeIniLocTxtFile(sub,hem)
%
% Create a text file of elec coordinates readable by NYU code. 
% output is voxel coordinates
%
% Author: David M. Groppe
% June, 2015

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
if strcmpi(hem(1),'l')
    postimpLocFname=sprintf('%s%sPostimpLocLeft.txt',elecReconPath,sub);
else
    postimpLocFname=sprintf('%s%sPostimpLocRight.txt',elecReconPath,sub);
end


%% space delimited from mgrid
[eCoords, elecLabels, elecRgb, elecPairs, elecPresent]=mgrid2matlab(sub,hem(1));
fprintf('Creating file: %s\n',postimpLocFname);
fid=fopen(postimpLocFname,'w');
for a=1:length(elecLabels)
    %elecLabels{a}=rmChar(elecLabels{a},'-'); % Hugh's electrode localization code no likey (fix later??)
    if ~isempty(findstr('grid',lower(elecLabels{a})))
        elecType='G';
    elseif ~isempty(findstr('depth',lower(elecLabels{a})))
        elecType='D';
        %elecLabels{a}=rmSubstring(elecLabels{a},'depth');
    else
        elecType='S';
    end
    fprintf(fid,'%s %f %f %f %s\n',elecLabels{a},eCoords(a,1), ...
        eCoords(a,2),eCoords(a,3),elecType);
    %     fprintf(fid,'%s %d %d %d %s\n',elecLabels{a},round(eCoords(a,1)), ...
    %         round(eCoords(a,2)),round(eCoords(a,3)),elecType);
end
fclose(fid);
