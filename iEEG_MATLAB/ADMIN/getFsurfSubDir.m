function fsurfSubDir = getFsurfSubDir()
%function fsurfSubDir = getFsurfSubDir()
%   Returns the Freesurfer Subject Directory as defined by the global
%   Matlab variable global_fs_dir or the shell variable SUBJECTS_DIR

global global_fs_dir;

if ~isempty(global_fs_dir)
    fsurfSubDir=global_fs_dir;
else
    if ispc,
        error('Hey mon, if you be using Windows you need to define the global variable "global_fs_dir" and put the path to your FreeSurfer subject folder in it.');
    else
        fsurfSubDir=getenv('SUBJECTS_DIR');
    end
end

end

