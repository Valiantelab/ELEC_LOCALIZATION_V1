%---------------------------------------------------
% FUNCTION: function loc_plot_mgh(electrodes,surface_path,v,colors,labels)
% INPUTS:   
%           electrodes = Nx3 matrix of RAS coordinates of electrodes
%           
%           surface_path = path to freesurfer output files (i.e. path to lh.pial and rh.pial)
%
%   	    v = 
%               'r' = view from right               or . . .
%               'l' = view from left                or . . .
%               'a' = view from anterior            or . . .
%               'p' = view from posterior           or . . .
%               's' = view from superior            or . . .
%               'i' = view from inferior            
%
%           colors = Nx3 vector of electrode colors
%               N = index number of electrode, each row of the Nx3 matrix is an rgb color value (e.g. [0 0 0] = black)
%
%           labels = Nx1 cell array of electrode labels (optional)
%     
%           
% OUTPUT:   Reconstructed cortical surface with electrode locations
%---------------------------------------------------
% last updated: 2010.03.31
%--------------------------------------------------- 

function loc_plot_mgh(electrodes,surface_path,v,colors,labels)

if nargin < 4
    colors = zeros(size(electrodes));
end

%read surfaces into matlab
[cortex.lh.vert, cortex.lh.tri] = read_surf([surface_path '/lh.pial']);
[cortex.rh.vert, cortex.rh.tri] = read_surf([surface_path '/rh.pial']);

%plot surfaces
h=gcf;
tripatch(cortex.rh,h,[.7 .7 .7]);
tripatch(cortex.lh,h,[.7 .7 .7]);
shading interp; lighting gouraud; material dull;

% adjust the view
switch v
    case 'r'
        l=light('Position',[1 0 0]);
        view(90,0)
    case 'l'
        l=light('Position',[-1 0 0]);
        view(270,0)
    case 'a'
        l=light('Position',[0 1 0]);  
        view(180,0)
    case 'p'
        l=light('Position',[0 -1 0]);        
        view(0,0)
    case 's'
        l=light('Position',[0 0 1]);    
        view(90,90)
    case 'i'
        l=light('Position',[0 0 -1]);
        view(270,270)
end
axis off, hold on


for j = 1:length(electrodes)
    plot3(electrodes(j,1),electrodes(j,2),electrodes(j,3),'o','Color',colors(j,:),'MarkerFaceColor', colors(j,:),'MarkerSize',10)
    hold all
    if nargin > 4
        h=text(electrodes(j,1),electrodes(j,2),electrodes(j,3),electrode_labels(j));
        set(h,'FontSize',6, 'Color', 'k');
    end
end

end
