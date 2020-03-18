function snp_insert_fig_attributes(doch,event)
%% 
%  
%  file:   snp_insert_fig_attributes.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2017.01.05. Thursday, 21:51:11
%

try
    
    fig = gcf;
    ax = gca;
    
    ret_copy = sprintf('''Position'', %s, ''Color'', [1 1 1]', pcz_num2str(fig.Position));
    
    ret_insert = sprintf('view(%s)\nset(gca,''Position'', %s)',...
        pcz_num2str(ax.View), pcz_num2str(ax.Position));
       
    disp(ret_copy)
    disp(ret_insert)

    Leg = findobj(fig, 'Type', 'Legend');
    if ~isempty(Leg)
        ret_insert = [ret_insert ...
            newline ...
            '% Leg = findobj(fig, ''Type'', ''Legend'');' newline ...
            sprintf('set(Leg,''Position'', %s, ''Box'', ''off'')', ...
                pcz_num2str(Leg.Position))
            ];
    end
    
    disp(ret_copy)
    disp(ret_insert)
    
    clipboard('copy', ret_copy)
    
    active = matlab.desktop.editor.getActive;
    editor = active.JavaEditor;
    editor.insertTextAtCaret(ret_insert);    
catch ex
    getReport(ex)
end

end