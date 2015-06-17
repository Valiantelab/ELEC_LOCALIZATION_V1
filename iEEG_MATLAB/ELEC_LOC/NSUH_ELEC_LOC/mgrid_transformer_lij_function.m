%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []= mgrid_transformer_lij_function(VOXEL_LIST,SURF_PATH,s,mgridfile,surftype)
% m.argyelan & sb 2010
% sb-25102011, added compute electrode distances
% sb-12142012, [pathstr, name, ext, versn]=fileparts(VOXEL_LIST) gave error
% too many outputs. changed to [pathstr, name, ext]

%VOXEL_LIST='/home/amiklos/work/VOYVODIC/electrode_localization/Miklos_fmr/Lobo_grids_freesurfer.VOX';
%SURF_PATH='	';
%s='r'; %side
%mgridfile='*.mgrid';
%surftype='pial'
[pathstr, name, ext]=fileparts(VOXEL_LIST);


if nargin < 5
    surftype='pial-outer-smoothed';
    pialoutfile=[pathstr '/' name '.PIALVOX'];
    undumpoutfile=[pathstr '/' name '.3dUndump.VOX'];
    picname=[pathstr '/' name];
elseif strcmpi(surftype,'pial')
    pialoutfile=[pathstr '/' name '.PIALVOX_nonSm'];
    undumpoutfile=[pathstr '/' name '.3dUndump_nonSm.VOX'];
    picname=[pathstr '/' name '_nonSm'];
elseif strcmpi(surftype,'white')
    pialoutfile=[pathstr '/' name '.PIALVOX_wm'];
    undumpoutfile=[pathstr '/' name '.3dUndump_wm.VOX'];
end
%keyboard
%VOX2RAS=[-1 0 0 128; 0 0 1 -128; 0 -1 0 128; 0 0 0 1];
VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];


%% make voxel list for later (used for *.pial.mgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VOX_bef=dlmread(VOXEL_LIST);
VOX_coor=VOX_bef(:,2:4);
RAS_coor=(VOX2RAS*[VOX_coor'; ones(1, size(VOX_coor,1))])';

%snap to dural (pial-outer-smoothed)
%RAS_aft=get_loc_snap_mgh(RAS_coor(:,1:3),SURF_PATH,s,surftype);
%keyboard
%snap to dural (pial-outer-smoothed)
if strcmpi(surftype,'pial')
    display('surftype=pial; first snap to smoothed pial, then pial.');
    RAS_aft_sm=get_loc_snap_mgh(RAS_coor(:,1:3),SURF_PATH,s,'pial-outer-smoothed');
    RAS_aft=get_loc_snap_mgh(RAS_aft_sm,SURF_PATH,s,surftype);
    get_vertex_nr=1; %for later to get txt file with vertex nr and surface distances
else
    RAS_aft=get_loc_snap_mgh(RAS_coor(:,1:3),SURF_PATH,s,surftype);
    get_vertex_nr=0;
end

RAS2VOX=inv(VOX2RAS);
VOX_aft=VOX_bef;
TEMP=(RAS2VOX*[RAS_aft'; ones(1, size(VOX_coor,1))])';
VOX_aft(:,2:4)=TEMP(:,1:3);

%% make list with electrode names (taken from mgrid) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ...and replace coordinates of depth electrodes with orig. mgrid coordinates
%(in previous step depth elecs where snapped to surface)

elec_names=mgrid_getNames(mgridfile);

% get index of depth electrodes
c=1; depth_ind=[];
for i=1:length(elec_names)
    if strfind(lower(elec_names{i}),'depth')
        depth_ind(c)=i;
        c=c+1;
    end
end
if ~isempty(depth_ind)
    VOX_aft(depth_ind,:)=VOX_bef(depth_ind,:);
end

%% write out pialvox and plot
%  dlmwrite([pathstr '/' name '.PIALVOX'],VOX_aft); %it writes out comma delimited
dlmwrite(pialoutfile,VOX_aft); %it writes out comma delimited
% red=[ones(1,size(RAS_coor,1));zeros(1,size(RAS_coor,1));zeros(1,size(RAS_coor,1))]';
if ~strcmpi(surftype,'white')
    if ~isempty(depth_ind); RAS_plot_noDepth=RAS_aft; RAS_plot_noDepth(depth_ind,:)=[];
    else RAS_plot_noDepth=RAS_aft;
    end
    red=[ones(1,size(RAS_plot_noDepth,1));zeros(1,size(RAS_plot_noDepth,1));zeros(1,size(RAS_plot_noDepth,1))]';
    loc_plot_mgh(RAS_plot_noDepth,SURF_PATH,s,red);
    print('-djpeg','-r300',picname);
end

%% make voxel-list readable by 3dUndump %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 columns: x / y / z / incrementing number (defines different value to each seed in 3dUndump)
VOX_aft_undump=zeros(size(VOX_aft));
VOX_aft_undump(:,1:2)=VOX_aft(:,2:3);
VOX_aft_undump(:,3)=255-VOX_aft(:,4);
VOX_aft_undump(:,4)=1:size(VOX_aft_undump,1);
dlmwrite(undumpoutfile,VOX_aft_undump,' '); %writes out space delimited

%% get vertex number and distance on surface
if get_vertex_nr==1; %run this only when snapping to pial surface
    display('get vertex numbers, replace depth elecs with zeros, write it to a txt file.');
    [vert_cart, vert_tri] = read_surf([SURF_PATH '/' s 'h.' surftype]);
    VERTEX_nr=zeros(1,size(RAS_aft,1));
    for i=1:size(RAS_aft,1)
        VERTEX_nr(i)=find(ismember(vert_cart,RAS_aft(i,:),'rows'))-1;
    end
    %replace depth electrodes with zeros
    if ~isempty(depth_ind)
        VERTEX_nr(depth_ind)=0;
    end
    dlmwrite([pathstr '/' name '.VertexNr'],VERTEX_nr,' ');
    
    
    %%% UNDER CONSTRUCTION %%%
    % get distance on pial surface and euclidian if depth
    % get fs subject name
    t=tokenize(SURF_PATH,'/');
    FS_subj=t{end-1};
    nruns=(length(VERTEX_nr)*length(VERTEX_nr))/2;
    
    use_euclid=1;
    use_sphere=0;
    use_fs=0;

    if use_euclid==1
        dist_euclid=zeros(length(VERTEX_nr));c=0;
        for i=1:length(VERTEX_nr)
            for ii=i+1:length(VERTEX_nr)
                c=c+1;
                coord1=VERTEX_nr(i); coord2=VERTEX_nr(ii);
                display([num2str(c) '/' num2str(nruns) ' compute euclidian distances for all electrode pairs.']);
                dist_euclid(i,ii)=sqrt( (RAS_aft(i,1)-RAS_aft(ii,1))^2 + (RAS_aft(i,2)-RAS_aft(ii,2))^2 + (RAS_aft(i,3)-RAS_aft(ii,3))^2 );
%                 dist_euclid(i,ii)=sqrt( (VOX_aft(i,1)-VOX_aft(ii,1))^2 + (VOX_aft(i,2)-VOX_aft(ii,2))^2 + (VOX_aft(i,3)-VOX_aft(ii,3))^2 );
            end
        end
        dlmwrite([pathstr '/' name '.distancesEuclid'],dist_euclid);
    end
    
    if use_sphere==1;
        dist_sphere=zeros(length(VERTEX_nr));c=0;
        for i=1:length(VERTEX_nr)
            for ii=i+1:length(VERTEX_nr)
                c=c+1;
                coord1=VERTEX_nr(i); coord2=VERTEX_nr(ii);
                if coord1==0 || coord2==0 %at least one elec is depth.
                    %compute euclidian distance.
                    display([num2str(c) '/' num2str(nruns) ' compute distance on sphere. but, at least one electrode is depth. get euclidian distance.']);
                    dist_sphere(i,ii)=sqrt( (RAS_aft(i,1)-RAS_aft(ii,1))^2 + (RAS_aft(i,2)-RAS_aft(ii,2))^2 + (RAS_aft(i,3)-RAS_aft(ii,3))^2 );
                else
                    %compute on sphere, DOUBLE CHECK
                    display([num2str(c) '/' num2str(nruns) ' compute distance on sphere.']);
                    [sphere_cart, sphere_tri] = read_surf([SURF_PATH '/' s 'h.sphere']);
                    dist_sphere(i,ii)=sphdist_from_cart(sphere_cart(coord1,:),sphere_cart(coord2,:));
                end
            end
        end
        dlmwrite([pathstr '/' name '.distancesSphere'],dist_sphere);
    end
    
    if use_fs==1; %(takes for freaking ever, 5-10h?)
        dist_surface=zeros(length(VERTEX_nr));c=0;
        for i=1:length(VERTEX_nr)
            for ii=i+1:length(VERTEX_nr)
                c=c+1;
                coord1=VERTEX_nr(i); coord2=VERTEX_nr(ii);
                if coord1==0 || coord2==0 %at least one elec is depth.
                    %compute euclidian distance. CHECK IF SURFACE AND THIS HAS SAME DIMENSION
                    display([num2str(c) '/' num2str(nruns) ' compute distance on surface. but, at least one electrode is depth. get euclidian distance.']);
                    dist_surface(i,ii)=sqrt( (RAS_aft(i,1)-RAS_aft(ii,1))^2 + (RAS_aft(i,2)-RAS_aft(ii,2))^2 + (RAS_aft(i,3)-RAS_aft(ii,3))^2 );
                    
                else %get distance on surface
                    display([num2str(c) '/' num2str(nruns) ' get distance on surface with mris_pmake (run in terminal)']);
                    fid = fopen(['/' pathstr '/run_getdist'], 'w');
                    t1=['mris_pmake --subject ' FS_subj  ' --hemi ' s 'h --surface0 pial --surface1 pial --curv0 curv --curv1 curv --mpmProg pathFind --mpmOverlay distance --mpmArgs startVertex:' num2str(coord1) ',endVertex:' num2str(coord2) '>' pathstr '/temp.txt' ];
                    fprintf(fid, '%s', t1);
                    fclose(fid);
                    system(['sh /' pathstr '/run_getdist']);
                    delete([pathstr '/run_getdist']);
		    %read out the distance value from temp.txt
                    fid=fopen([pathstr '/temp.txt']); 
                    tline = fgetl(fid);
                    tline = fgetl(fid);
                    fclose(fid);
                    x=tokenize(tline,'[]');
                    dist_surface(i,ii)=str2double(x{2});
                end
            end
        end
        dlmwrite([pathstr '/' name '.distancesFS'],dist_surface);
    end
end
