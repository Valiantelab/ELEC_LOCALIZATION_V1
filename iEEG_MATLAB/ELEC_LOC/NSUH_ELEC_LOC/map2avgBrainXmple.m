% Right hem subs
groupAvgCoords=[];
groupLabels=[];
patientLabels=[];
subs={'TWH001','TWH002'};
ct=0;
for a=1:length(subs),
    [avgCoords, avg_vids, sub_vids, labels]=pvox2AvgBrain(subs{a},'r',[]);
    groupAvgCoords=[groupAvgCoords; avgCoords];
    groupLabels=[groupLabels; labels];
    for b=1:size(avgCoords,1),
        ct=ct+1;
        patientLabels{ct}=subs{a};
    end
end