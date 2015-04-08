function clickText(h,txt)
%function clickText(h,txt)
%
% Adds code to a figure object so that when you click on it, a text box
% appears with the desired text in it. When you click on the text box, it
% disappears.
%
% Inputs
%  h   - figure object handle (vector or singleton)
%  txt - string of cell array of strings containing object labels
%
% Author: David Groppe
% Mehtalab, 2012

if iscell(txt)
    if length(h)~=length(txt)
       error('To the number of elements of h and txt are different.'); 
    end
    for a=1:length(h),
        set(h(a),'userdata',txt{a});
        bdfcn=['Cp = get(gca,''CurrentPoint''); ' ...
            'Xp=Cp(1,1);', ...
            'Yp=Cp(1,2);', ...
            'Zp=Cp(1,3);', ...
            'dat=get(gcbo,''userdata'');', ...
            'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));', ...
            'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'');'];
        set(h,'buttondownfcn',bdfcn);
    end
else
    set(h,'userdata',txt);
    bdfcn=['Cp = get(gca,''CurrentPoint''); ' ...
        'Xp=Cp(1,1);', ...
        'Yp=Cp(1,2);', ...
        'Zp=Cp(1,3);', ...
        'dat=get(gcbo,''userdata'');', ...
        'ht=text(Xp,Yp,Zp,sprintf(''%s'',dat));', ...
        'set(ht,''backgroundcolor'',''w'',''horizontalalignment'',''center'',''verticalalignment'',''middle'',''buttondownfcn'',''delete(gcbo);'');'];
    set(h,'buttondownfcn',bdfcn);
end