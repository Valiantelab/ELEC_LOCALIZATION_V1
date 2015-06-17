% plotElecPial() - Function for plotting Freesurfer surfaces
%                            with or without colored maps or electrodes.
%
% Usage:
%  >> cfg_out=plotElecPial(sub_name,cfg);
%
% Required Input:
%   sub_name - Name of the subject's freesurfer directory (full path not
%              needed)
%
% Optional Inputs:
%   cfg variable with the following possible fields:
%
%    Electrode Options:
%     eleccoord            -If 'n', no electrodes will be rendered in the
%                           figure.  If 'y', electrode coordinates will be
%                           taken from *.pialvox file in patient's
%                           FreeSurfer folder.  Alternatively, you
%                           can pass a 2D matrix of coordinates
%                           instead. {default: 'y'}
%     elecsize             -Size of electrode markers (disks or spheres).
%                           This also determines thickness of lines connecting
%                           electrodes (if any) and electrode labels (if
%                           shown). {default=8};
%     elecshape            -'marker' or 'sphere': The shape used to
%                           represent electrodes. {default: 'marker'}
%     eleccolors           -2D matrix of colors to fill electrodes
%                           (rows=electrodes, columns=RGB values) or a vector 
%                           of values that will be automatically converted 
%                           into a colorscale.
%                           {default: all electrodes filled with black}
%     edgeblack            -If 'y', electrodes will all have a black
%                           border. Otherwise, border will be same color as
%                           marker. This argument has no effect if
%                           electrodes are represented as spheres. {default:
%                           'y'}
%     elecnames            -Cell array of the names of the electrodes to
%                           which the rows of eleccolors corresponds.
%                           Electrodes not included in elecnames will be
%                           colored black. Note, THESE SHOULD BE MGRID
%                           NAMES. {default: not used}
%     linewidth            -Thickness of line connecting pairs of
%                           electrodes. {default: elecsize/3}
%     badchans             -Cell array of the names of bad electrodes that
%                           should be plotted black (for when other
%                           electrodes are in color). {default: not used}
%     ignorechans          -Cell array of the names of electrodes that
%                           will not be shown. {default: not used}
%     ignoredepthelec    -'y' or 'n': If 'y', depth electrodes will not
%                           be shown. Depth electrodes have the string
%                           "depth" in their mgrid file names. {default: 'y'}
%     onlyshow             -Cell array of the names of the only electrodes
%                           you want to be shown. If an electrode is listed
%                           both in onlyshow and ignorechans, it will not
%                           be shown. {default: not used}
%     pullout              -Factor via which to project electrodes out from
%                           the center of view. Helpful for when electrodes
%                           sink into the cortical surface. {default: 1}
%     showlabels           -'y' on 'n': If 'y', the name of the first and
%                           last electrode of each strip and each corner of
%                           the 64 chan grid will be shown next to electrode.
%                           {default: 'y'}
%     pairs                -A nx3, 4 or 5 cell array specifying n pairs of
%                           electrodes to be connected with lines.
%                           The first two columns indicate which electrodes
%                           are in the pair.
%                           The third column is a 3 element vector indicating
%                           the RGB color of the line that will be drawn
%                           to join the pair of electrodes.
%                           The fourth, optional, column is the text that
%                           will appear when the line joining the electrodes
%                           is clicked on.
%                           The fifth, optional, column is the connection
%                           strength for each pair. The maximum value will
%                           correspond to linewidth; others will be a ratio
%                           of that value.
%                           {default: not used}
%     plotcbar             -'y' or 'n': Plot colorbar next to brain. {default:
%                           'y' if funcfname eleccolors argument specified, 
%                           'n' otherwise}
%     colorscale           -'absmax','minmax', 'justpos', 'justneg',
%                           or numeric vector [minval maxval].  The limits
%                           that define the electrode data color scale.
%                           {default: 'absmax'}.
%     units                -A string or []. The title of the colorbar
%                           (e.g., 'z-score'). If empty, not title is drawn.
%                           {default: []}
%
%    Surface Options:
%     surftype             -'pial' or 'inflated': Type of Freesurfer surface
%                           to plot. If inflated, gyri and sulci are marked
%                           dark and light grey. {default='pial'};
%     overlay_parcellation -If 'DK', Freesurfer Desikan-Killiany (36 area)
%                           cortical parcellation is plotted on the surface.
%                           If 'D', Freesurfer Destrieux (76 area)
%                           parcellation is used. {default: not used}
%     colors_parcellation  -Optional N-by-3 numeric array with RGB indexes
%                           (0:255) for each of the ROIs in the
%                           overlay_parcellation. The colors for each ROI
%                           need to be in the exact same order as that of
%                           the default Freesurfer color table (you can get
%                           that by using the Freesurfer-MATLAB function
%                           [averts,albl,actbl=read_annotation(fname.annot);
%                           the ROIs are listed in actbl.struct_names and
%                           their numeric labels (found in albl) in
%                           actbl.table).
%                           {default: not used; the default Freesurfer
%                           colors for the parcellation are used instead}
%     opaqueness           -[1 to 0] The "alpha" level of the pial surface.
%                           0 means the surface is completely transparent.
%                           1 means that it is completely opaque.
%                           {default: 1}
%     view                 -Angle and lighting with which to view brain.
%                           This also defines which hemisphere to plot.
%                           Options include:
%                             'omni' - 6 views of each hemisphere. When you
%                             use this option with the funcfname option,
%                             funcfname needs to be a cell array of two
%                             filenames. The first specifies the left hem
%                             values and the second the right.
%                             'lomni' - 6 views of left hemisphere
%                             'l' - Left hem lateral
%                             'lm' - Left hem medial
%                             'lo' - Left hem occipital
%                             'lf' - Left hem frontal
%                             'lim' - Left hem inferior-medial
%                             'li' - Left hem inferior
%                             'ls' - Left hem superior
%                             'lsv' - Left hem superior and vertically
%                                     aligned
%                             'liv' - Left hem inferior and vertically
%                                     aligned
%                           Replace 'l' with 'r' to get these views of the
%                           right hemisphere. Alternatively, you can define
%                           everything yourself like so:
%                                   brain_view.light
%                                   brain_view.hem
%                                   light: [1 0 0]
%                                   eyes: [45 0]
%                                   hem: 'r'
%
%
%    Other Options:
%     axis                 -Handle of axis in which to make plot.
%                           {default: new axis created}
%     figid                -Handle of figure in which to make plot.
%                           {default: new figure created}
%     clearfig             -'y' or 'n'; If 'y', the figured is cleared
%                           before the plot is created. {default: 'y??'}
%     backgroundcolor      -Standard Matlab color argument (e.g., 'k' or
%                           [.1 .5 .3]). The axis background color.
%                           {default: not used}
%     rotate3d             -'y' or 'n'; If 'y', clicking on the axis will
%                           allow user to manually rotate it.  If 'n',
%                           clicking on the electrode will reveal the
%                           electrode's name.  Clicking again on the
%                           electrode name will remove it. {default: ??}
%     title                -Title to place on figure. {default: 'y'}
%     fsurfsubdir          -The path to the FreeSurfer subject directory.
%                           Necessary if running MATLAB on Windows.
%                           {default: taken from shell}
%     clearglobal          -If 'n', some plotting information is left in
%                           global memory. Useful for speeding *omni plots.
%                           {default: 'y'}
%     verblevel            - An integer specifying the amount of information you want
%                           this function to provide about what it is doing during runtime.
%                            Options are:
%                             0 - quiet, only show errors, warnings, and external function reports
%                             1 - stuff anyone should probably know
%                             2 - stuff you should know the first time you start working
%                                 with a data set {default value}
%                             3 - stuff that might help you debug (show all
%                                 reports)
%
% Example:
% % Plot electrodes on brain with Desikan-Killiany cortical parcellation
% cfg=[];
% cfg.view='r';
% cfg.figid=1;
% cfg.overlay_parcellation='DK';
% cfg.rotate3d='n';
% cfg.showlabels='y';
% cfg.title=[];
% cfg_out=plotElecPial('TWH014',cfg);
%
% % Plot depths with a semi-transparent pial surface
% cfg=[];
% cfg.view='ri';
% cfg.ignoredepthelec='n';
% cfg.opaqueness=.5;
% cfg.title=[];
% cfg_out=plotElecPial('TWH014',cfg);
%
% % Plot electrodes as spheres, color coded to reflect correlation value
% elecnames=cell(6,1);
% for a=1:6,
%    elecnames{a}=sprintf('RMF%d',a); 
% end
% cfg=[];
% cfg.view='romni';
% cfg.figid=1;
% cfg.elecshape='sphere';
% cfg.eleccolors=rand(6,1);
% cfg.colorscale='minmax';
% cfg.showlabels='n';
% cfg.units='r';
% cfg.elecnames=elecnames;
% cfg.title='TWH014: Pie Man Correlations';
% cfg_out=plotElecPial('TWH014',cfg);
%
%  Authors:
%   David M. Groppe, Stephan Bickel, Pierre Megevand, Andrew Dykstra
%   Laboratory for Multimodal Human Brain Mapping
%   Feinstein Institute for Medical Research
%   Manhasset, New York
%


% HISTORY:
% adykstra/ck 08-2010; get_loc_snap_mgh
% sb 05-2011: do not save graph option
% sb 06-2011: added config structure as input
% dg 03-2011: massively re-written
% dg 5/15/2011: ESM bars between electrodes now also affected by "pullout"
%  option
% dg 1/3/2013: Updated comments and got rid of underscores in most arguments.
%  Also slightly modified SB's modifications to funcfname argument so that
%  you can pass a vector of values instead of a filename.
% dg 5/13/13: Comments to "pairs" option added (and now pairs can have four
%  columns). click_text changed to clickText3D so that text doesn't
%  disappear into brain.
% dg 9/9/13 changed nominal colors from Destrieux palette to
%  distinguishable_colors.m 
% dg 1/2/14 omni, lomni, and romni views added; cfg.units option added
% mm 24/02/14: pairs can know have 5 columns to specify linewidth for each
%  pairs
% pm x/xx/14 debug passing 2D numeric array for eleccoord
% pm 7/18/14 added optional colors_parcellation input
% pm 7/24/14 debug passing 2D numeric array for eleccoord AND overlay_parcellation
% pm 7/28/14 debug passing 2D numeric array for eleccoord AND inflated pial surface
% pm 8/12/14 debug ?omni view
% dg 4/15 eleccolors can now be a vector of values from which a colorscale
% is automatically derived
% dg 6/15 now expects electrode coordinates and names in Yang method format

%% TO DO
% do auto colormap for electrodes (both jet and rb)?
% Size of electrodes should be relative to brain size?
% Make eleccolors and colorbar work for bipolar lines too

function cfg_out=plotElecPial(subj,cfg)

%% Parse parameters
if ~isfield(cfg, 'elecsize'),       elecsize = 8;          else  elecsize = cfg.elecsize;      end
if ~isfield(cfg, 'snap2surf'),      snap2surf = 0;         else  snap2surf = cfg.snap2surf;      end
if ~isfield(cfg, 'surftype'),       surftype = 'pial';     else  surftype = cfg.surftype;     end
if ~isfield(cfg, 'eleccoord'),      eleccoord= 'yes';      else  eleccoord = cfg.eleccoord;       end
if ~isfield(cfg, 'eleccolors'),     eleccolors= [];        else  eleccolors = cfg.eleccolors;        end
if ~isfield(cfg, 'plotcbar'),       plotcbar=[];          else plotcbar=cfg.plotcbar; end
if ~isfield(cfg, 'units'),          units=[];            else units=cfg.units; end
if ~isfield(cfg, 'colorscale'),     colorscale='absmax';   else colorscale=cfg.colorscale; end
if ~isfield(cfg, 'pullout'),        pullout= 1;            else  pullout = cfg.pullout; end
if ~isfield(cfg, 'view'),           brain_view= 'l';       else  brain_view = cfg.view; end
if ~isfield(cfg, 'axis'),           h_ax=[];               else  h_ax=cfg.axis; end
if ~isfield(cfg, 'overlay_parcellation'), overlay_parcellation=0;  else  overlay_parcellation=cfg.overlay_parcellation; end
if ~isfield(cfg, 'colors_parcellation'),  colors_parcellation= []; else  colors_parcellation= cfg.colors_parcellation;  end
if ~isfield(cfg, 'figid'),         h_fig=[];              else  h_fig=cfg.figid; end
if ~isfield(cfg, 'clearfig'),       clearfig=1;            else  clearfig=cfg.clearfig; end
if ~isfield(cfg, 'title'),          surf_title='default';  else surf_title=cfg.title; end
if ~isfield(cfg, 'elecnames'),      color_elecnames=[];    else color_elecnames=cfg.elecnames; end
if ~isfield(cfg, 'badchans'),       badchans=[];           else badchans=cfg.badchans; end
if ~isfield(cfg, 'ignorechans'),    ignorechans=[];        else ignorechans=cfg.ignorechans; end
if ~isfield(cfg, 'rotate3d'),       rotate3d_active='yes'; else rotate3d_active=cfg.rotate3d; end
if ~isfield(cfg, 'fsurfsubdir'),  fs_dir=[];             else fs_dir=cfg.fsurfsubdir; end
if ~isfield(cfg, 'verblevel'),      verblevel=2;           else verblevel=cfg.verblevel; end
if ~isfield(cfg, 'backgroundcolor'), backgroundcolor=[]; else backgroundcolor=cfg.backgroundcolor; end
if ~isfield(cfg, 'pairs'), electrode_pairs=[];  else electrode_pairs=cfg.pairs; end
if ~isfield(cfg, 'linewidth'), linewidth=[];  else linewidth=cfg.linewidth; end
if ~isfield(cfg, 'onlyshow'), onlyshow=[];  else onlyshow=cfg.onlyshow; end
if ~isfield(cfg, 'ignoredepthelec'), ignoredepthelec='y'; else ignoredepthelec=cfg.ignoredepthelec; end
if ~isfield(cfg, 'edgeblack'),      edgeblack='y';         else edgeblack=cfg.edgeblack; end
if ~isfield(cfg, 'elecshape'), electrodeshape='marker';  else electrodeshape=cfg.elecshape; end
if ~isfield(cfg, 'showlabels'), showlabels=0;  else showlabels=universalYes(cfg.showlabels); end
if ~isfield(cfg, 'opaqueness'),      opaqueness=1;          else opaqueness=cfg.opaqueness; end
if ~isfield(cfg, 'clearglobal'),    clearglobal=1;          else clearglobal=cfg.clearglobal; end

checkCfg(cfg,'plotElecPial.m');
cmapName=[];

if strcmpi(cfg.view,'omni')
    cfg_out=plotPialOmni(subj,cfg);
    return;
elseif strcmpi(cfg.view,'lomni') || strcmpi(cfg.view,'romni')
    cfg_out=plotPialHemi(subj,cfg);
    return;
end


global global_fs_dir;
if isempty(fs_dir)
    if ~isempty(global_fs_dir)
        fs_dir=global_fs_dir;
    else
        if ispc,
            error('Hey mon, if you be using Windows you need to be specifying the fsurfsubdir input argument.');
        else
            fs_dir=getenv('SUBJECTS_DIR');
        end
    end
end
% Folder with surface files
surfacefolder=fullfile(fs_dir,subj,'surf');
if ~isempty(surfacefolder) && (surfacefolder(end)~='/')
    surfacefolder=[surfacefolder '/'];
end

% Get side of brain to show
if ischar(brain_view),
    if findstr(brain_view,'l')
        side='l';
    else
        side='r';
    end
else
    side=lower(brain_view.hem);
    if ~strcmpi(side,'l') && ~strcmpi(side,'r')
        error('cfg.brain_view.hem needs to be ''l'' or ''r''.');
    end
end

if isempty(eleccolors)
    if isempty(plotcbar)
        plotcbar='n';
    end;
else
    if isempty(plotcbar)
        plotcbar='y';
    end;
end


%% %%%%%%%%%% Start Main Function %%%%%%%
VerbReport('**** PLOTTING CORTICAL SURFACE WITH "plotElecPial.m" ****', ...
    2,verblevel);


%% MAKE FIGURE
if ~isempty(h_fig),
    figure(h_fig);
else
    h_fig=figure;
end
if universalYes(clearfig),
    clf;
end
if ~isempty(backgroundcolor)
    set(h_fig,'color',backgroundcolor);
end
if ~isempty(h_ax),
    axes(h_ax);
else
    h_ax=gca;
end


%% If plotting on inflated surface, load curvature values so that sulci and
% gyri can be seen
if strcmpi(surftype,'inflated')
    if side == 'r'
        [curv, fnum] = read_curv([surfacefolder 'rh.curv']);
    else
        [curv, fnum] = read_curv([surfacefolder 'lh.curv']);
    end
    curv_map=zeros(length(curv),3);
    pcurv_ids=find(curv>=0);
    curv_map(pcurv_ids,:)=repmat([1 1 1]*.3,length(pcurv_ids),1);
    ncurv_ids=find(curv<0);
    curv_map(ncurv_ids,:)=repmat([1 1 1]*.7,length(ncurv_ids),1);
end
%global map cbar_min cbar_max;
global cbar_min cbar_max;
if strcmp(surftype,'inflated');
    map=curv_map;
else
    map=[.7 .7 .7]; %make it all grey
end


%% READ SURFACE
global cort %speeds up omni a tiny bit
if isempty(cort)
    if side == 'r'
        [cort.vert cort.tri]=read_surf([surfacefolder 'rh.' surftype]);
    else
        [cort.vert cort.tri]=read_surf([surfacefolder 'lh.' surftype]);
    end
    if min(min(cort.tri))<1
        cort.tri=cort.tri+1; %sometimes this is needed sometimes not. no comprendo. DG
    end
end


%% PLOT SURFACE
n_map_vert=size(map,1);
tripatchDG(cort,h_fig,map); %this plots the brain
if universalYes(rotate3d_active),
    rotate3d on;
else
    rotate3d off;
end


%% If specified, overlay cortical parcellation
if overlay_parcellation,
    labelfolder=fullfile(fs_dir,subj,'label');
    if ~isempty(labelfolder) && (labelfolder(end)~='/')
        labelfolder=[labelfolder '/'];
    end
    if strcmpi(overlay_parcellation,'DK')
        annot_fname=[labelfolder side 'h.aparc.annot']; %Desikan-Killiany 36 area atlas
        [averts,albl,actbl]=read_annotation(annot_fname);
    elseif strcmpi(overlay_parcellation,'D')
        annot_fname=[labelfolder side 'h.aparc.a2009s.annot']; %Destrieux 76 area atlas
        [averts,albl,actbl]=read_annotation(annot_fname);
        actbl.table(43,1:3)=255*[1 1 1]*.7; %make medial wall the same shade of grey as functional plots
    else
        error('overlay_parcellation argument needs to take a value of ''D'' or ''DK''.');
    end
    if ~isempty(colors_parcellation)
        if size(colors_parcellation,1)~=size(actbl.table,1)
            error('plotElecPial:colors_parcellation_size1','colors_parcellation argument needs to have \n the same number of rows as the number of ROIs \n in the parcellation. For %s, %d.',overlay_parcellation,size(actbl.table,1));
        end
        if size(colors_parcellation,2)~=3
            error('plotElecPial:colors_parcellation_size2','colors_parcellation must be an N-by-3 array.');
        end
        actbl.table(:,1:3)=colors_parcellation;
    end
    clear averts;
    [~,loc_table]=ismember(albl,actbl.table(:,5));
    loc_table(loc_table==0)=1; % for the case where the label for the vertex says 0
    fvcdat=actbl.table(loc_table,1:3)./255;
    clear loc_table;
    tsurf_h=trisurf(cort.tri, cort.vert(:, 1), cort.vert(:, 2), cort.vert(:, 3),...
        'FaceVertexCData', fvcdat,'FaceColor', 'interp','FaceAlpha',1);
    if ~strcmpi(eleccoord,'no') && ~strcmpi(eleccoord,'n')
        cfg_elec2parc=[];
        cfg_elec2parc.fsurfsubdir=fs_dir;
        if isnumeric(eleccoord)
            cfg_elec2parc.eleccoord=eleccoord;
            cfg_elec2parc.elecnames=cfg.elecnames;
        end
        elec_assign=[];
    end
else
    elec_assign=[];
end

%% Set Lighting & View
shading interp; lighting gouraud; material dull; axis off, hold on
if ischar(brain_view)
    switch brain_view
        case 'r'
            l=light('Position',[1 0 0]);
            view(90,0)
        case 'rm'
            l=light('Position',[-1 0 0]);
            view(270,0)
        case 'rim'
            l=light('Position',[-1 0 0]);
            view(270,-45)
        case 'ri'
            l=light('Position',[0 0 -1]);
            view(90,-90)
            %case 'ru'
            %    l=light('Position',[1 0 0]);
            %    view(270,-90)
        case 'ro'
            l=light('Position',[0 -1 0]);
            view(0,0)
        case 'lo'
            l=light('Position',[0 -1 0]);
            view(0,0)
        case 'rf'
            l=light('Position',[0 1 0]);
            view(180,0)
        case 'lf'
            l=light('Position',[0 1 0]);
            view(180,0)
        case 'rs'
            l=light('Position',[0 0 1]);
            view(90,90);
        case 'rsv' %superior & vertically aligned
            l=light('Position',[0 0 1]);
            view(0,90);
        case 'l'
            l=light('Position',[-1 0 0]);
            view(270,0);
        case 'lm'
            l=light('Position',[1 0 0]);
            view(90,0);
        case 'li'
            l=light('Position',[0 0 -1]);
            view(90,-90);
        case 'lim'
            l=light('Position',[-1 0 0]);
            view(270,-45);
            %case 'lu'
            %    l=light('Position',[-1 0 0]);
            %    view(270,-90);
        case 'ls'
            l=light('Position',[0 0 1]);
            view(270,90);
        case 'lsv' %superior & vertically aligned
            l=light('Position',[0 0 1]);
            view(0,90);
        case 'liv' %inferior & vertically aligned
            l=light('Position',[0 0 -1]);
            view(0,-90);
        case 'riv' %inferior & vertically aligned
            l=light('Position',[0 0 -1]);
            view(0,-90);
    end
    clear l
else
    light('Position',brain_view.light);
    view(brain_view.eyes)
end
alpha(opaqueness);


%% PLOT ELECTRODES (optional)
elec_sphere=0; %default
if strcmpi(eleccoord,'no') || strcmpi(eleccoord,'n')
    display('...not plotting electrodes');
else
    %VOX2RAS=[-1 0 0 128; 0 0 -1 128; 0 -1 0 128; 0 0 0 1];
    if isnumeric(eleccoord) % PM moved that section up here
        if size(eleccoord,2)==3
            display('...Electrode input is matrix with coordinates.');
            RAS_coor=eleccoord;
            elec_names=cfg.elecnames; % PM added this line; uses cfg.elecnames on purpose!
            if strcmpi(surftype,'inflated')
                cfg_pvox2inf=[];
                cfg_pvox2inf.fsurfsubdir=fs_dir;
                cfg_pvox2inf.eleccoord=RAS_coor;
                cfg_pvox2inf.elecnames=elec_names;
                RAS_coor=pvox2InfBrain(subj,side,cfg_pvox2inf);
            end
        else
            error('...Electrode input is numeric but doesn''t have 3 coordinates');
        end
    else % electrode coordinates and names to be read from .DURAL and .electrodeNames files
        if universalYes(eleccoord)
            display('...Overlaying electrodes. Taking *.DURAL from elec_recon folder. Use cfg.eleccord=''n''; if not wanted.');
            if strcmpi(side,'L')
                f=dir([fs_dir '/' subj '/elec_recon/*left.DURAL']);
                if isempty(f)
                    error('You need to specify electrode coordinates or have a file that matches: %s\n', ...
                        [fs_dir '/' subj '/elec_recon/*left.DURAL']);
                end
            else
                f=dir([fs_dir '/' subj '/elec_recon/*right.DURAL']);
                if isempty(f)
                    error('You need to specify electrode coordinates or have a file that matches: %s\n', ...
                        [fs_dir '/' subj '/elec_recon/*right.DURAL']);
                end
            end
            if length(f)>1,
                error('Too many possible DURAL files.  I do not know which to use.');
            end
            eleccoord=[fs_dir '/' subj '/elec_recon/' f(1).name];
        end
        
        %if exist(eleccoord,'file') && findstr(eleccoord,'VOX')
        if exist(eleccoord,'file')
            if strcmpi(surftype,'inflated')
                cfg_pvox2inf=[];
                cfg_pvox2inf.fsurfsubdir=fs_dir;
                RAS_coor=pvox2InfBrain(subj,side,cfg_pvox2inf);
            else
                tempCsv=csv2Cell(eleccoord,' ',2);
                RAS_coor=zeros(size(tempCsv,1),1);
                for csvLoopA=1:size(tempCsv,1),
                    for csvLoopB=1:3,
                        RAS_coor(csvLoopA,csvLoopB)=str2num(tempCsv{csvLoopA,csvLoopB});
                    end
                end
            end
            [path, nme, ext]=fileparts(eleccoord);
            if strcmpi(side,'L')
                f=dir([fs_dir '/' subj '/elec_recon/*left.electrodeNames']);
                if isempty(f)
                    error('You need to specify electrode names or have a file that matches: %s\n', ...
                        [fs_dir '/' subj '/elec_recon/*left.electrodeNames']);
                end
            else
                f=dir([fs_dir '/' subj '/elec_recon/*right.electrodeNames']);
                if isempty(f)
                    error('You need to specify electrode names or have a file that matches: %s\n', ...
                        [fs_dir '/' subj '/elec_recon/*right.electrodeNames']);
                end
            end
            if length(f)>1,
                error('Too many possible electrodeNames files.  I do not know which to use.');
            end
            
            enames_filename=[fs_dir '/' subj '/elec_recon/' f(1).name];
            fprintf('Attempting to read electrode names from file %s\n', ...
                enames_filename);
            elec_names=csv2Cell(enames_filename,' ',1);
            elec_names=format_elec_names(elec_names);
        end
    end
    
    % exclude depth electrodes
    if universalYes(ignoredepthelec),
        c=1; depth_ind=[];
        for i=1:length(elec_names)
            if ~isempty(strfind(lower(elec_names{i}),'depth')) || strcmpi(elec_names{i}(1),'d')
                depth_ind(c)=i;
                c=c+1;
            end
        end
        if ~isempty(depth_ind),
            fprintf('%d depth electrodes removed.\n',length(depth_ind));
            RAS_coor(depth_ind,:)=[];
            temp_elecnames=cell(1,1);
            c=0;
            for i=1:length(elec_names),
                if ~ismember(i,depth_ind),
                    c=c+1;
                    temp_elecnames{c}=elec_names{i};
                end
            end
            elec_names=temp_elecnames;
            clear temp_elecnames;
        end
    end
    
    % Pull electrodes out from the brain towards the viewer
    nRAS=size(RAS_coor,1);
    if pullout,
        fprintf('...pulling out electrodes by factor %f. cfg.pullout=0 if not wanted.\n',pullout);
        v=axis;
        campos=get(gca,'cameraposition');
        %camtarg=get(gca,'cameratarget'); ?? Should we check that this be set to 0?
        err=repmat(campos,nRAS,1)-RAS_coor;
        nrmd=err./repmat(sqrt(sum(err.^2,2)),1,3);
        RAS_coor=RAS_coor+nrmd*pullout;
    end
    
    % Plot lines joining electrodes
    if ~isempty(electrode_pairs),
        electrode_pairs(:,1)=format_elec_names(electrode_pairs(:,1));
        electrode_pairs(:,2)=format_elec_names(electrode_pairs(:,2));
        if isempty(linewidth)
            linewidth=elecsize/3;
        end
        if size(electrode_pairs,2)<=4
            electrode_pairs(:,5) ={linewidth};
        elseif size(electrode_pairs,2)>4
            % normalize linewidth
            electrode_pairs(:,5) = cellfun(@rdivide,electrode_pairs(:,5),num2cell(repmat(max([electrode_pairs{:,5}]), [size(electrode_pairs,1) 1])),'UniformOutput',false);
            electrode_pairs(:,5) = cellfun(@times,electrode_pairs(:,5),num2cell(repmat(linewidth, [size(electrode_pairs,1) 1])),'UniformOutput',false);
        end
        
        n_pairs=size(electrode_pairs,1);
        pair_ids=[0 0];
        for a=1:n_pairs,
            for b=1:2,
                [got_it, pair_ids(b)]=ismember(lower(electrode_pairs{a,b}),lower(elec_names));
                if ~got_it
                    error('Channel %s is in electrode pairs but not in pialVox electrode names.',electrode_pairs{a,b});
                end
            end
            hl=plot3([RAS_coor(pair_ids(1),1) RAS_coor(pair_ids(2),1)], ...
                [RAS_coor(pair_ids(1),2) RAS_coor(pair_ids(2),2)], ...
                [RAS_coor(pair_ids(1),3) RAS_coor(pair_ids(2),3)],'-');
            if isnumeric(electrode_pairs{a,3})
                set(hl,'color',electrode_pairs{a,3},'linewidth',linewidth);
            else
                set(hl,'color',str2num(electrode_pairs{a,3}),'linewidth',linewidth);
            end
            
            if size(electrode_pairs,2)>3
                clickText3D(hl,[electrode_pairs{a,1} '-' electrode_pairs{a,2} ': ' electrode_pairs{a,4}],2);
            else
                clickText3D(hl,[electrode_pairs{a,1} '-' electrode_pairs{a,2}],2);
            end
        end
    end
    
    % make electrodes black if no input given
    if isempty(eleccolors)
        eleccolors = zeros(size(RAS_coor));
    elseif isvector(eleccolors)
        if isnumeric(cfg.colorscale)
            type='minmax';
            cbar_min=cfg.colorscale(1);
            cbar_max=cfg.colorscale(2);
        else
            type=cfg.colorscale;
        end
        if verLessThan('matlab','8.0.1')
            cmapName='jet';
        else
            cmapName='parula';
        end
        if isempty(cbar_min),
            % make electrode colormap
            [eleccolors, elecLimits, cmapName]=vals2Colormap(eleccolors,type,cmapName);
            cbar_min=elecLimits(1); 
            cbar_max=elecLimits(2);
        else
            [eleccolors, elecLimits, cmapName]=vals2Colormap(eleccolors,type,cmapName,[cbar_min cbar_max]);
        end
    else
       % DG ?? to do, grab cbar_min and cbar_max from optional inputs when matrix of rgb values passed cbar_min=min(
       
    end
    if ~isempty(color_elecnames),
        color_elecnames=format_elec_names(color_elecnames);
        n_color_electrodes=length(color_elecnames);
        used_color_electrodes=zeros(1,n_color_electrodes);
    end
    if isempty(onlyshow),
        %if user didn't specify a subset of electrodes to show, attempt to
        %show all of them
        onlyshow=elec_names;
    else
        onlyshow=format_elec_names(onlyshow);
    end
    if ~isempty(ignorechans),
        onlyshow=setdiff(onlyshow,ignorechans);
    end
    
    % Prepare variables if electrodes are to be drawn as spheres
    if strcmpi(electrodeshape,'sphere')
        elec_sphere=1;
        [sphX, sphY, sphZ]=sphere(elecsize*6);
        Zdim=size(sphZ);
        scale_sph=2;
        sphX=sphX*scale_sph;
        sphY=sphY*scale_sph;
        sphZ=sphZ*scale_sph;
        sph_colors=zeros(nRAS,3);
        sph_ct=0;
    else
        elec_sphere=0;
    end
    for j = 1:nRAS
        if ismember(lower(elec_names{j}),lower(onlyshow)),
            if ~isempty(color_elecnames)
                [have_color id]=ismember(lower(elec_names{j}),lower(color_elecnames));
                if ~have_color || (~isempty(badchans) && ismember(lower(elec_names{j}),lower(badchans)))
                    % We don't have data for this electrode or have
                    % been instructed to ignore it; plot it black
                    if elec_sphere
                        sph_ct=sph_ct+1;
                        h_elec(sph_ct)=surf(sphX+RAS_coor(j,1),sphY+RAS_coor(j,2),sphZ+RAS_coor(j,3),zeros(Zdim));
                        %h_sph(sph_ct)=surf(sphX+RAS_coor(j,1),sphY+RAS_coor(j,2),sphZ+RAS_coor(j,3),zeros(Zdim));
                        sph_colors(sph_ct,:)=[1 1 1]*.01;
                        %colormap(gca,[1 1 1]*.01);
                        %shading interp; lighting gouraud; material dull;
                        %set(h_sph(sph_ct),'facecolor',[1 1 1]*.01);
                    else
                        h_elec=plot3(RAS_coor(j,1),RAS_coor(j,2),RAS_coor(j,3),'o','Color','k','MarkerFaceColor','k','MarkerSize',elecsize);
                    end
                    if showlabels,
                        add_name(RAS_coor(j,:),elec_names{j},elec_names,elecsize,[1 1 1])
                    end
                else
                    % Color the electrode to represent data
                    if elec_sphere
                        sph_ct=sph_ct+1;
                        h_elec(sph_ct)=surf(sphX+RAS_coor(j,1),sphY+RAS_coor(j,2),sphZ+RAS_coor(j,3),zeros(Zdim));
                        sph_colors(sph_ct,:)=eleccolors(id,:);
                        used_color_electrodes(id)=1;
                        %colormap(gca,eleccolors(id,:));
                        %shading interp; lighting gouraud; material dull;
                        %set(h_sph(j),'facecolor',eleccolors(id,:));
                    else
                        if universalYes(edgeblack)
                            markeredgecolor=[0 0 0];
                        else
                            markeredgecolor=eleccolors(id,:);
                        end
                        h_elec=plot3(RAS_coor(j,1),RAS_coor(j,2),RAS_coor(j,3),'o', ...
                            'Color',eleccolors(id,:),'MarkerFaceColor', eleccolors(id,:),'MarkerSize',elecsize, ...
                            'MarkerEdgeColor',markeredgecolor,'linewidth',2);
                        used_color_electrodes(id)=1;
                    end
                    if showlabels,
                        add_name(RAS_coor(j,:),elec_names{j},elec_names,elecsize,eleccolors(id,:))
                    end
                end
            else
                if elec_sphere
                    sph_ct=sph_ct+1;
                    h_elec(sph_ct)=surf(sphX+RAS_coor(j,1),sphY+RAS_coor(j,2),sphZ+RAS_coor(j,3),zeros(Zdim));
                    %colormap(gca,eleccolors(j,:))
                    %set(h_sph(j),'facecolor',eleccolors(j,:));
                    %shading interp; lighting gouraud; material dull;
                    sph_colors(sph_ct,:)=eleccolors(j,:);
                else
                    h_elec=plot3(RAS_coor(j,1),RAS_coor(j,2),RAS_coor(j,3),'o','Color',eleccolors(j,:),'MarkerFaceColor', eleccolors(j,:),'MarkerSize',elecsize);
                    if showlabels,
                        add_name(RAS_coor(j,:),elec_names{j},elec_names,elecsize,eleccolors(j,:))
                    end
                end
            end
            hold all
            if ~universalYes(rotate3d_active),
                if isempty(elec_assign)
                    set(h_elec,'userdata',elec_names{j});
                else
                    set(h_elec,'userdata',[elec_names{j} ' ' elec_assign{j,2}]);
                end
                % This click_text code should put the text out towards the
                % viewer (so it doesn't get stuck in the brain)
                pop_fact=5; %this might be too far for lateral surfaces
                bdfcn=['Cp = get(gca,''CurrentPoint''); ' ...
                    'Cp=Cp(1,1:3);', ...
                    'v=axis;', ...
                    'campos=get(gca,''cameraposition'');', ...
                    'df=Cp-campos;', ...
                    'nrmd=df/sqrt(sum(df.^2));', ...
                    sprintf('Cp=Cp-%d*nrmd;',pop_fact), ...
                    'dat=get(gcbo,''userdata'');', ...
                    'ht=text(Cp(1),Cp(2),Cp(3),sprintf(''%s'',dat));', ...
                    'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'');'];
                set(h_elec,'buttondownfcn',bdfcn);
            end
        end
        %NOTE:
        % x dimension is lateral/medial  (+=lateral)
        % y dimension is ant/posterior (+=anterior)
        % z dimension is superior/inferior (+=superior)
    end
    if verblevel>1,
        if ~isempty(color_elecnames),
            not_found=find(used_color_electrodes==0);
            if not_found
                warning('Number of colored electrodes NOT FOUND: %d\n',length(not_found));
                for dg=not_found,
                    fprintf('%s\n',color_elecnames{dg});
                end
            end
        end
    end
    if elec_sphere,
        shading interp; lighting gouraud; material dull;
        %for some reason the shading command resets of the colors of all the
        %spheres, thus the need for this silly loop.  There is surely a more
        %elegant way to deal with this.
        for a=1:sph_ct,
            set(h_elec(a),'facecolor',sph_colors(a,:));
        end
    end
end


%% Add Title
if strcmpi(surf_title,'default'),
    if ischar(brain_view)
        surf_title=[subj '; ' brain_view '; '];
    else
        surf_title=[subj '; ' brain_view.hem '; '];
    end
    title(surf_title,'fontsize',20);
elseif ~isempty(surf_title)
    title(surf_title,'fontsize',20);
end



%% Colorbar
h_cbar=[]; % needs to be declared for cfg_out even if colorbar not drawn
if universalYes(plotcbar)
    if isempty(cbar_min) || isempty(cbar_max)
        fprintf('cbar_min or cbar_max are empty. Cannot draw colorbar.\n');
    else
        nColors=64; % # of colors in jet colormap
        %colormap(cmapName);
        h_cbar=cbarDG('vert',[1:nColors],[cbar_min cbar_max],5,cmapName);
        
        if isequal(get(h_fig,'color'),[0 0 0]);
            %If background of figure is black, make colorbar text white
            set(h_cbar,'xcolor','w'); % fix so that box isn't white? ??
            set(h_cbar,'ycolor','w');
        end
        if ~isempty(units)
            ht=title(units); % ?? change size?
        end
    end
end

%% COLLECT CONFIG OUTPUT
% Update this? DG ??
cfg_out.subject=subj;
cfg_out.view=brain_view;
cfg_out.elecsize=elecsize;
cfg_out.surftype=surftype;
cfg_out.h_cbar=h_cbar;
cfg_out.h_brain=h_ax;
cfg_out.cmapName=cmapName;
if exist('cfg','var'), cfg_out.cfg=cfg; end
if exist('RAS_coor','var'), cfg_out.electrode_coords=RAS_coor; end
if exist('cbar_min','var'), cfg_out.cbar_limits=[cbar_min cbar_max]; end

if universalYes(clearglobal)
    clear global cbar_min cbar_max cort;
end

%%%% END OF MAIN FUNCTION %%%%


%% HELPER FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function elec_names=format_elec_names(elec_names)
% Removes underscores from electrode names and replaces 'Grid' with 'G'
n_chan=length(elec_names);
for a=1:n_chan,
    % Remove underscores
    n_char=length(elec_names{a});
    temp_str=repmat('a',1,n_char);
    ct=0;
    for b=1:n_char,
        if elec_names{a}(b)~='_'
            ct=ct+1;
            temp_str(ct)=elec_names{a}(b);
        end
    end
    elec_names{a}=temp_str(1:ct);
    if length(elec_names{a})>4
        if strcmpi(elec_names{a}(1:4),'Grid')
            elec_names{a}=['G' elec_names{a}(5:end)];
        end
    end
end


%% subfunction ADD_NAME
function add_name(xyz,label,all_labels,markersize,rgb)
% right now rgb argument is ignored, fix in the future so that electrode
% names stand out from background? ??

if keyElec(label,all_labels),
    h_t=text(xyz(1),xyz(2),xyz(3),label);
    set(h_t,'color','k','fontweight','bold','fontsize',markersize+4);
    %     if isequal(rgb,[0 0 0]),
    %         set(h_t,'color','w','fontweight','bold','fontsize',markersize+2);
    %     else
    %         set(h_t,'color','k','fontweight','bold','fontsize',markersize+2);
    %     end
end


%% subfunction ODD_ELEC
function yesno=odd_elec(elec_name)
%function yesno=odd_elec(elec_name) ?? used?

num_ids=find((elec_name>=48).*(elec_name<=57));
elec_num=str2num(elec_name(num_ids));
if mod(elec_num,2)
    yesno=1;
else
    yesno=0;
end


%% subfunction plotPialOmni
function sub_cfg_out=plotPialOmni(subj,cfg)

if ~isfield(cfg, 'figid'),         h_fig=[];            else  h_fig=cfg.figid; end
if ~isfield(cfg, 'zthresh'),       zthresh=[];          else  zthresh = cfg.zthresh; end
if ~isfield(cfg, 'figid'),         h_fig=[];              else  h_fig=cfg.figid; end
if ~isfield(cfg, 'fsurfsubdir'),   fs_dir=[];             else fs_dir=cfg.fsurfsubdir; end
if ~isfield(cfg, 'eleccoord'),      eleccoord= 'yes';      else  eleccoord = cfg.eleccoord;       end
if ~isfield(cfg, 'elecsize'),       elecsize = 8;          else  elecsize = cfg.elecsize;      end
if ~isfield(cfg, 'eleccolors'),     eleccolors= [];        else  eleccolors = cfg.eleccolors;        end
if ~isfield(cfg, 'colorscale'),     colorscale=[];   else colorscale=cfg.colorscale; end
if ~isfield(cfg, 'units'),     units=[];   else units=cfg.units; end
if ~isfield(cfg, 'showlabels'),         showlabels='y';            else  showlabels=cfg.showlabels; end
if ~isfield(cfg, 'plotcbar'),     plotcbar=[];   else plotcbar=cfg.plotcbar; end

if ~isempty(eleccolors),
    if isempty(plotcbar)
        plotcbar='y';
    end
end

global global_fs_dir;
if isempty(fs_dir)
    if ~isempty(global_fs_dir)
        fs_dir=global_fs_dir;
    else
        if ispc,
            error('Hey mon, if you be using Windows you need to be specifying the fsurfsubdir input argument.');
        else
            fs_dir=getenv('SUBJECTS_DIR');
        end
    end
end

clear global cbar_min cbar_max cort;

if isempty(h_fig),
    h_fig=figure; clf;
else
    figure(h_fig); clf;
end
set(h_fig,'MenuBar','none','position',[100 190 1000 600],'paperpositionmode','auto');

%% Figure out which hemisphere has electrodes
f_left=dir([fs_dir '/' subj '/elec_recon/*left.DURAL']);
left_coverage=~isempty(f_left);
f_right=dir([fs_dir '/' subj '/elec_recon/*right.DURAL']);
right_coverage=~isempty(f_right);

%%
used_limits=[];
for h=1:2,
    for v=1:6, %will run 1-6
        ax_loc=[0 0 0 0];
        if h==1,
            bview='l';
            %h_ax(1)=axes('position',[.03 .01 .45 .95]);
        else
            bview='r';
        end
        if v==2,
            bview=[bview 'm'];
        elseif v==3,
            bview=[bview 'f'];
        elseif v==4,
            bview=[bview 'o'];
        elseif v==5,
            bview=[bview 'sv'];
        elseif v==6,
            bview=[bview 'iv'];
        end
        switch (h-1)*6+v,
            case 1 %LH lateral
                ax_loc=[0 .67 .4 .3];
            case 2 %LH medial
                ax_loc=[0 .34 .4 .3];
            case 3 %LH frontal
                ax_loc=[.155 .02 .2 .3];
            case 4 %LH occiptal
                ax_loc=[1-.15-.2 .02 .2 .3];
            case 5 %LH superior
                ax_loc=[1-.455-.2 .55 .2 .4];
            case 6 %LH inferior
                ax_loc=[1-.455-.2 .14 .2 .4];
                
            case 7 %%%%% RH lateral
                ax_loc=[1-.4 .67 .4 .3];
            case 8 %RH medial
                ax_loc=[1-.4 .34 .4 .3];
            case 9 %RH frontal
                ax_loc=[.045 .02 .2 .3];
            case 10 %RH occipital
                ax_loc=[1-.035-.2 .02 .2 .3];
            case 11 %RH superior
                ax_loc=[.455 .55 .2 .4];
            case 12 %RH inferior
                ax_loc=[.455 .14 .2 .4];
        end
        fprintf('CODE %d, VIEW %s\n',(h-1)*6+v,bview);
        h_ax=axes('position',ax_loc);
        sub_cfg=cfg;
        sub_cfg.view=bview;
        
        sub_cfg.title=[];
        if bview(1)=='l'
            if ~left_coverage || (ischar(eleccoord) && universalNo(eleccoord)), % PM added ischar for eleccoord
                sub_cfg.eleccoord='n';
            else
                sub_cfg.eleccoord='y';
                if ~isfield('elecsize',sub_cfg)
                    sub_cfg.elecsize=6;
                end
                sub_cfg.showlabels=showlabels;
            end
        else
            if ~right_coverage || (ischar(eleccoord) && universalNo(eleccoord)), % PM added ischar for eleccoord
                sub_cfg.eleccoord='n';
            else
                sub_cfg.eleccoord='y';
                if ~isfield('elecsize',sub_cfg)
                    sub_cfg.elecsize=6;
                end
                sub_cfg.showlabels=showlabels;
            end
        end
        
        sub_cfg.figid=h_fig;
        sub_cfg.axis=h_ax;
        
        sub_cfg.clearfig='n';
        sub_cfg.rotate3d='n';
        sub_cfg.plotcbar='n';
        if v==6
            sub_cfg.clearglobal=1; %last view for this hem, clear map from global memory
        else
            sub_cfg.clearglobal=0;
        end
        sub_cfg_out=plotElecPial(subj,sub_cfg);
        
        if isempty(used_limits)
            if isfield(sub_cfg_out,'cbar_limits')
                used_limits=sub_cfg_out.cbar_limits;
            end
        else
            if used_limits(2)<sub_cfg_out.cbar_limits(2)
                used_limits(2)=sub_cfg_out.cbar_limits(2);
            end
            if used_limits(1)>sub_cfg_out.cbar_limits(1)
                used_limits(1)=sub_cfg_out.cbar_limits(1);
            end
        end
    end
end

if universalYes(plotcbar),
    % Colorbar
    h_cbar=axes('position',[.4 .06 .2 .03]);
    map=colormap;
    n_colors=size(map,1);
    n_tick=5;
    cbarDG(h_cbar,[1:n_colors],used_limits,n_tick,sub_cfg_out.cmapName);
    if ~isempty(units)
        ht=title(units);
        set(ht,'fontsize',12);
    end
    ticklabels=cell(1,n_tick);
    ticks=linspace(used_limits(1),used_limits(2),n_tick);
    for a=1:n_tick,
        ticklabels{a}=num2str(ticks(a));
    end
    set(h_cbar,'yticklabel',ticklabels);
end

if isfield(cfg,'title')
    if ~isempty(cfg.title)
        % Overall Fig Title
        ht=textsc(cfg.title,'title');
        set(ht,'fontweight','bold','fontsize',20,'position',[0.5 0.975]);
    end
end
drawnow;

%% subfunction plotPialHemi
function sub_cfg_out=plotPialHemi(subj,cfg)

if ~isfield(cfg, 'usemask'),       usemask=[];          else usemask=cfg.usemask; end
if ~isfield(cfg, 'figid'),         h_fig=[];            else  h_fig=cfg.figid; end
if ~isfield(cfg, 'fsurfsubdir'),   fs_dir=[];             else fs_dir=cfg.fsurfsubdir; end
if ~isfield(cfg, 'eleccoord'),      eleccoord= 'yes';      else  eleccoord = cfg.eleccoord;       end
if ~isfield(cfg, 'elecsize'),       elecsize = 8;          else  elecsize = cfg.elecsize;      end
if ~isfield(cfg, 'eleccolors'),     eleccolors= [];        else  eleccolors = cfg.eleccolors;        end
if ~isfield(cfg, 'colorscale'),     colorscale=[];   else colorscale=cfg.colorscale; end
if ~isfield(cfg, 'units'),     units=[];   else units=cfg.units; end
if ~isfield(cfg, 'plotcbar'),     plotcbar=[];   else plotcbar=cfg.plotcbar; end

if ~isempty(eleccolors),
    if isempty(plotcbar)
        plotcbar='y';
    end
end

global global_fs_dir;
if isempty(fs_dir)
    if ~isempty(global_fs_dir)
        fs_dir=global_fs_dir;
    else
        if ispc,
            error('Hey mon, if you be using Windows you need to be specifying the fsurfsubdir input argument.');
        else
            fs_dir=getenv('SUBJECTS_DIR');
        end
    end
end


clear global cbar_min cbar_max cort;

hem=cfg.view(1);

if isempty(h_fig),
    h_fig=figure; clf;
else
    figure(h_fig); clf;
end
set(h_fig,'MenuBar','none','position',[100 190 800 500],'paperpositionmode','auto');

%% Viewpoints
used_limits=[];
for v=1:6, %will run 1-6
    ax_loc=[0 0 0 0];
    bview=hem;
    if v==2,
        bview=[bview 'm'];
    elseif v==3,
        if bview(1)=='r',
            bview=[bview 'f'];
        else
            bview=[bview 'o'];
        end
    elseif v==4,
        if bview(1)=='r'
            bview=[bview 'o'];
        else
            bview=[bview 'f'];
        end
    elseif v==5,
        bview=[bview 's'];
    elseif v==6,
        bview=[bview 'i'];
    end
    switch v,
        case 1 % lateral
            ax_loc=[-.05 .52 .55 .45];
        case 2 % medial
            ax_loc=[-.05 .05 .55 .45];
        case 3 % frontal
            ax_loc=[.29 .55 .55 .41];
        case 4 % occiptal
            ax_loc=[.56 .55 .44 .41];
        case 5 % superior
            ax_loc=[.41 .05 .54 .23];
        case 6 % inferior
            ax_loc=[.41 .31 .54 .23];
    end
    fprintf('CODE %d, VIEW %s\n',v,bview);
    h_ax=axes('position',ax_loc);
    
    sub_cfg=cfg;
    sub_cfg.view=bview;
    sub_cfg.title=[];
    sub_cfg.figid=h_fig;
    sub_cfg.axis=h_ax;
    sub_cfg.clearfig='n';
    sub_cfg.rotate3d='n';
    sub_cfg.plotcbar='n';
    if v==6
        sub_cfg.clearglobal=1; %last view, clear plotElecPial from global memory
    else
        sub_cfg.clearglobal=0;
    end
    sub_cfg_out=plotElecPial(subj,sub_cfg);
    
    if isempty(used_limits)
        if isfield(sub_cfg_out,'cbar_limits')
            used_limits=sub_cfg_out.cbar_limits;
        end
    else
        if used_limits(2)<sub_cfg_out.cbar_limits(2)
            used_limits(2)=sub_cfg_out.cbar_limits(2);
        end
        if used_limits(1)>sub_cfg_out.cbar_limits(1)
            used_limits(1)=sub_cfg_out.cbar_limits(1);
        end
    end
end

%% Colorbar
if universalYes(plotcbar),
    h_cbar=axes('position',[.90 .1 .03 .8]);
    map=colormap;
    n_colors=size(map,1);
    n_tick=5;
    cbarDG(h_cbar,1:n_colors,used_limits,n_tick,sub_cfg_out.cmapName);
    if ~isempty(units)
        ht=title(units);
        set(ht,'fontsize',12);
    end
    ticklabels=cell(1,n_tick);
    ticks=linspace(used_limits(1),used_limits(2),n_tick);
    for a=1:n_tick,
        ticklabels{a}=num2str(round(10*ticks(a))/10);
    end
    set(h_cbar,'yticklabel',ticklabels);
end


if isfield(cfg,'title')
    if ~isempty(cfg.title)
        % Overall Fig Title
        ht=textsc(cfg.title,'title');
        set(ht,'fontweight','bold','fontsize',20,'position',[0.5 0.975]);
    end
end
drawnow;

%subfunction VerbReport
function VerbReport(report,verbtag,VERBLEVEL)
% VerbReport() - Outputs messages if they exceed a desired level of importance.  
%
% Usage:
%  >> VerbReport(report,verbtag,VERBLEVEL)
%
% Inputs:
%   report    = a string that is some message to the user
%   verbtag   = an intger specifiying the importance of report
%   VERBLEVEL = an integer specifiying a threshold of importance
%               for displaying reports. If verbtag is less than VERBLEVEL, the
%               report will be displayed..
%
% Author: Tom Urbach
% Kutaslab

if nargin<3
  tmpVERBLEVEL = 3;
elseif isempty(VERBLEVEL)
  tmpVERBLEVEL = 3;
else
  tmpVERBLEVEL = VERBLEVEL;
end;

if verbtag <= tmpVERBLEVEL
	if ischar(report)
		fprintf('%s\n',report);
	else
		fprintf('%d\n',report);
	end;
end;
