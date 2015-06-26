% A simple script for testing the mapping to the freesufr avg brain

%% Collect coordinates on average brain
groupAvgCoords=[];
groupLabels=[];
subs={'TWH011','TWH014'};
%subs={'TWH011'};
for a=1:length(subs),
    [avgCoords, avg_vids, sub_vids, labels]=pvox2AvgBrain(subs{a},'l',[]);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; labels];
end

%% Plot electrodes
cfg=[];
cfg.view='l';
cfg.figid=10;
cfg.eleccoord=groupAvgCoords;
cfg.elecnames=groupLabels;
cfg.rotate3d='n';
cfg_out=plotElecPial('fsaverage',cfg);

%%
cfg=[];
cfg.view='l';
cfg.figid=10;
cfg.eleccoord=avgCoords;
cfg.elecnames=labels;
cfg.surftype='inflated';
cfg_out=plotElecPial('fsaverage',cfg);
%cfg_out=plotPial('fsaverage',cfg);