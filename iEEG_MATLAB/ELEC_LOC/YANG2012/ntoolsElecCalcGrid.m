function [out_cell, elec_stats, info_cell]=ntoolsElecCalcGrid(ini_cell,subjectpath,scale,radius,nRow,nCol)
%function [out_cell, elec_stats, info_cell]=ntoolsElecCalcGrid(ini_cell,subjectpath,scale,radius,nRow,nCol)
%
% calculate the grid electrodes with initial locations, ini_cell is a cell
% array in which first column is the [elec_name number], the rest 3
% columns are the x y z coordinates
%
% output: elec is a structure contains the grid names and their locations,

if isempty(ini_cell)
    disp('No grid electrode.');
    elec = []; elec_stats = []; info_cell = [];
    return;
end

ini_pos=zeros(4,1); % electrode #s
ini_loc=cell2mat(ini_cell(:,2:4)); % electrode initial RAS coordinates
for a=1:4,
    tempId=find(ini_cell{a,1}=='_');
    nameStem=ini_cell{a,1}(1:(tempId-1));
    ini_pos(a)=str2num(ini_cell{a,1}(tempId+1:end));
end

% determine the hemisphere that grid locates
if ini_loc(:,1)>0
    sph = 'rh';
elseif ini_loc(:,1)<0
    sph = 'lh';
else
    error('wrong initial positions for grid');
end

s=[nRow nCol]; % Grid dimensions

% project the grid on the outer-brain surface
[elec_proj, info_cell]= ntools_elec_projection(ini_loc,ini_pos,s(1),s(2),sph,subjectpath,scale,radius);
elec_stats = elec_proj(1,1:4);

% I think info_cell is just elec_proj at each iteration and elec_proj is
% the first iteration

% clear the unnecessary data
clear elec_temp elec_num ini_pos ini_loc;
%clear elec_temp elec_num elec_proj ini_pos ini_loc;
fprintf('Grid stats: itr#: %f mean: %f std:%f\n\n',cell2mat(elec_stats(1:3)));


%% get the data from the struct
nElec=nRow*nCol;
out_cell=cell(nElec,5);
for a=1:nElec,
    out_cell{a,1}=[nameStem '_' num2str(a)]; %Name fix this ??
    out_cell{a,2}=elec_proj{5}(a,1); %R
    out_cell{a,3}=elec_proj{5}(a,2); %A
    out_cell{a,4}=elec_proj{5}(a,3); %S
    out_cell{a,5}='G';
end

% l = 1;
% for j = 1:size(elec_temp2.(char(name{i})),1)
%     name_num(l) = cellstr(sprintf('%s%.2d',char(name{i}),j));
%     elec_pos(l,:) = elec_temp2.(char(name{i}))(j,:);
%     l = l+1;
% end
% 
% elec = cell(size(name_num,2),5);
% elec(:,1) = upper(name_num)';
% elec(:,2:4) = num2cell(elec_pos);
% 
% elec(:,5) = repmat({'G'},[l-1 1]);

fprintf('Done \n\n');