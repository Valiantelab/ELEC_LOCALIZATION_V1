function interpStripElec(fSub,hem)
%UNTITLED Summary of this function goes here
%   
% This function interpolates the location of a strip electrode by:
% 1) Using its neighbors to define a line
% 2) Finding the point on the line at the known 2D distance between the
% contacts
% 3) Projecting that point to the nearest smoothed dural surface vertex
%
% Information is taken from and output to the elec_recon subfolder of the
% patient's FreeSurfer folder.
%
% Author: David Groppe
% Honeylab, University of Toronto
% June, 2015

% Future work:
% -Edit the localization log to indicate that this function was run
% -If ever useful, make code able to combine interpolations from multiple
% neigbhors for mid-strip or grid electrodes.

global global_fs_dir;
if ~isempty(global_fs_dir)
    fsDir=global_fs_dir;
else
    fsDir=getenv('SUBJECTS_DIR');
end
if isempty(fsDir)
    error('Freesurfer subject directory not defined. Use global variable "global_fs_dir" to do this.');
end

% Print out possible electrodes
if strcmpi(hem,'r') || strcmpi(hem,'rh')
    hemLong='right';
elseif strcmpi(hem,'l') || strcmpi(hem,'lh')
   hemLong='left'; 
else
    error('Invalid value for "hem" argument.');
end
if fsDir(end)=='/'
    erPath=[fsDir fSub '/elec_recon/'];
    surfPath=[fsDir fSub '/surf/'];
else
    erPath=[fsDir '/' fSub '/elec_recon/'];
    surfPath=[fsDir '/' fSub '/surf/'];
end
inFname=[erPath fSub '_' hemLong '.electrodeNames'];
fprintf('Loading file %s\n',inFname);
eNames=csv2Cell(inFname,' ',1);


%% Choose electrode to interpolate
nElec=size(eNames,1);
elecStem=cell(nElec,1);
elecNum=zeros(nElec,1);
for a=1:nElec,
    id=find(eNames{a}=='_');
    elecStem{a}=eNames{a}(1:(id-1));
    elecNum(a)=str2num(eNames{a}(id+1:end));
end
uniStem=unique(elecStem);
nUni=length(uniStem);
fprintf('Electrode Strip/Grids:\n');
for a=1:nUni,
    %uniStem{a}=rmSubstring(uniStem{a},'depth');
    fprintf('%d: %s\n',a,uniStem{a});
end
stemId=NaN;
while isnan(stemId)
    stemId=input('Enter the number of strip/grid with mis-localized contact:');
    if stemId<1 || stemId>nUni,
       disp('Invalid response.');
       stemId=NaN;
    end
end
psblChans=findstrInCell(uniStem{stemId},elecStem);
nPsblChans=length(psblChans);
if nPsblChans<3,
   error('This function cannot fix a strip with less than 3 contacts.');
end
fprintf('Possible channels:\n');
for a=1:nPsblChans,
   fprintf('%d: %s%d\n',a,elecStem{psblChans(a)},elecNum(psblChans(a))); 
end
ctctId=NaN;
while isnan(ctctId)
    ctctId=input('Enter the number of strip/grid with mis-localized contact:');
    if ctctId<1 || ctctId>nPsblChans,
        disp('Invalid response.');
        ctctId=NaN;
    end
end
ctctNum=elecNum(psblChans(ctctId));
srtdPsblChans=sort(elecNum(psblChans));
if ctctNum==srtdPsblChans(nPsblChans)
    nborA=ctctNum-1;
    nborB=ctctNum-2;
else
    nborA=ctctNum+1;
    nborB=ctctNum+2;
end
nborAlong=[uniStem{stemId} '_' num2str(nborA)];
nborBlong=[uniStem{stemId} '_' num2str(nborB)];
fixLong=[uniStem{stemId} '_' num2str(ctctNum)];
fprintf('Using dural position of channels %s and %s to infer location of %s\n', ...
    nborAlong,nborBlong,fixLong);

%% Define line via neighbors
% import coordinates
inFname=[erPath fSub '_' hemLong '.DURAL'];
fprintf('Loading dural RAS coordinates from file %s\n',inFname);
eRAS=csv2Cell(inFname,' ',2);
if length(eRAS)~=length(eNames)
   error('# of electrode names ~= # of electrode coordinates'); 
end
ras=zeros(nElec,3);
for a=1:nElec,
    for b=1:3,
       ras(a,b)=str2num(eRAS{a,b}); 
    end
end
fixId=findstrInCell(fixLong,eNames,1);
nborIdA=findstrInCell(nborAlong,eNames,1);
nborIdB=findstrInCell(nborBlong,eNames,1);

% Project to a point down a line
neoRas=ras;
neoRas(fixId,:)=ras(nborIdA,:)+(ras(nborIdA,:)-ras(nborIdB,:));

%% Snap to nearest smoothed dural surface
if strcmpi(hemLong,'right')
    inSurf=[surfPath 'rh.pial-outer-smoothed'];
else
    inSurf=[surfPath 'lh.pial-outer-smoothed'];
end
duralSurf=fs_read_surf(inSurf);
if size(duralSurf.vertices,1)==3,
    % This might vary with version of freesurfer code
   duralSurf.vertices=duralSurf.vertices'; 
end
vertId = dsearchn(duralSurf.vertices,neoRas(fixId,:));
neoRas(fixId,:)=duralSurf.vertices(vertId,:);


%%
% Load CT coorindates for comparison
inFname=[erPath fSub '_' hemLong '.CT'];
fprintf('Loading CT coordinates from file %s\n',inFname);
ctRAS=csv2Cell(inFname,' ',2);
if length(ctRAS)~=length(eNames)
   error('# of electrode names ~= # of electrode coordinates'); 
end
rasCt=zeros(nElec,3);
for a=1:nElec,
    for b=1:3,
       rasCt(a,b)=str2num(ctRAS{a,b}); 
    end
end

%%
plotCtVsDural(neoRas,rasCt,eNames,fSub,hem,[fsDir '/' fSub],0,0);
% higlight changed electrode
plot3(neoRas(fixId,1),neoRas(fixId,2),neoRas(fixId,3),'g*');
hT=title(sprintf('Red=Postcorrect, Blue=Precorrect, Green=%s',rmChar(fixLong,'_')));

fprintf('Is interpolated location reasonable?\n');
resp=binaryQuery();
if universalYes(resp)
    %output new coordinates to file
    
    % RAS COORDINATES
    % Dural
    fnameDuralRAS = [erPath fSub '_' hemLong '.DURAL'];
    fprintf('Saving dural RAS electrode locations to: %s\n',fnameDuralRAS);
    fid=fopen(fnameDuralRAS,'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'R A S\n');
    for a=1:nElec,
        fprintf(fid,'%f %f %f\n',neoRas(a,1),neoRas(a,2),neoRas(a,3));
    end
    fclose(fid);
    
    % Pial
    %% snap to pial
    if strcmpi(hemLong,'right')
        inSurf=[surfPath 'rh.pial'];
    else
        inSurf=[surfPath 'lh.pial'];
    end
    pialSurf=fs_read_surf(inSurf);
    if size(pialSurf.vertices,1)==3,
        % This might vary with version of freesurfer code
        pialSurf.vertices=pialSurf.vertices';
    end
    pialRAS=zeros(nElec,3);
    for a=1:nElec,
        vertId = dsearchn(pialSurf.vertices,neoRas(a,:));
        pialRAS(a,:)=pialSurf.vertices(vertId,:);
    end
    fnamePialRAS = [erPath fSub '_' hemLong '.PIAL'];
    fprintf('Saving pial RAS electrode locations to: %s\n',fnamePialRAS);
    fid=fopen(fnamePialRAS,'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'R A S\n');
    for a=1:nElec,
        fprintf(fid,'%f %f %f\n',pialRAS(a,1),pialRAS(a,2),pialRAS(a,3));
    end
    fclose(fid);
    
    %% Re-output VOX coordinates
    % VOX COORDINATES
    VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];
    RAS2VOX=inv(VOX2RAS);
    duralVOX=(RAS2VOX*[neoRas'; ones(1, nElec)])';
    fnameDuralVOX = [erPath fSub '_' hemLong '.DURALVOX'];
    fprintf('Saving dural VOX electrode locations to: %s\n',fnameDuralVOX);
    fid=fopen(fnameDuralVOX,'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'X Y Z\n');
    for a=1:nElec,
        fprintf(fid,'%f %f %f\n',duralVOX(a,1),duralVOX(a,2),duralVOX(a,3));
    end
    fclose(fid);
    
    pialVOX=(RAS2VOX*[pialRAS'; ones(1, nElec)])';
    fnamePialVOX = [erPath fSub '_' hemLong '.PIALVOX'];
    fprintf('Saving pial VOX electrode locations to: %s\n',fnamePialVOX);
    fid=fopen(fnamePialVOX,'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'X Y Z\n');
    for a=1:nElec,
        fprintf(fid,'%f %f %f\n',pialVOX(a,1),pialVOX(a,2),pialVOX(a,3));
    end
    fclose(fid);
    
    %% Replot distance figs
    close;
    plotCtVsDural(neoRas,rasCt,eNames,fSub,hem,[fsDir '/' fSub],1,1); 
    
    % add fix to diary
    fprintf('Add note to 00README.txt file to indicate that this command was run on %s\n',datestr(now));
else
    fprintf('Exiting. Electrode coordinates NOT modified.\n');
end


end

