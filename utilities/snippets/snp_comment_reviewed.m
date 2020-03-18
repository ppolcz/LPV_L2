function returnval = snp_comment_reviewed(doch,event)
%%
%
%  file:   snp_comment_reviewed.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.31. Sunday, 12:17:19
%  Reviewed on 2016.02.17. Wednesday, 10:45:50
%  Reviewed on 2016.03.13. Sunday, 12:51:26
%  Reviewed on 2017. August 25. [date be informative]
%  Reviewed on 2017. September 24. [just return (.mlx)]
%  Modified on 2018. March 03. (Modified instead of Reviewed)
%

ret = sprintf('%%  Major review on %s (%s)\n', pcz_fancyDate('informative'), version('-release'));

try
    active = matlab.desktop.editor.getActive;
    editor = active.JavaEditor;

    linenr = editor.getLineNumber+1;
    lbegin = editor.lineAndColumnToPosition(linenr,0);
    editor.setCaretPosition(lbegin);
    editor.insertTextAtCaret(ret);    
catch ex
    getReport(ex)
    returnval = ret;
    
    ret = sprintf('Major review on %s (%s)', pcz_fancyDate('informative'), version('-release'));
    clipboard('copy',ret);
end

end

