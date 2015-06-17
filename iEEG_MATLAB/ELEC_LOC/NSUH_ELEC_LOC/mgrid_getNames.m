function [electrode_names]=mgrid_getNames(mgridfile)

% read mgrid file. get name of electrodes.
% sb-05062010
% sb-12142012, fileparts gave error "too many outputs" deleted "vers"

[pathstr, name, ext]=fileparts(mgridfile);
file_new=[pathstr '/' name '.electrodeNames'];

if exist(file_new,'file')
    delete(file_new);
end

fid =fopen(mgridfile);
fid2 =fopen(file_new,'a');
counter=1;
last_elec='';


while  1
    tline = fgetl(fid)
    if ~ischar(tline) break, end
    if strfind(tline,'#Description')
        elec_name = fgetl(fid); 
        if ~strcmp(elec_name,'patient')
            while ~strcmp(tline,'#Dimensions')
                tline=fgetl(fid);
            end
            dim=fgetl(fid);
            dimn=str2num(dim);
            %get electrode number (eg for a 4x5 grid in sequence in mgrid is 
            % written 5/10/15/20/4/9...) 
            elec_numbers=ones(dimn(1),dimn(2));
            c=0;
            for i=1:dimn(2)
                elec_numbers(:,i)=(dimn(2)-c):dimn(2):(dimn(2)*dimn(1)-c);
                c=c+1;
            end
            
            numel=dimn(1)*dimn(2);
            for i=1:numel
               out=[elec_name '_' num2str(elec_numbers(i))];
               strfind(out,'depth')
               electrode_names{counter}=out;
               counter=counter+1;
               fprintf(fid2,'%s\n',out);
            end
        end
    end
end

fclose(fid)
fclose(fid2)

