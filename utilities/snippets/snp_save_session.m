function [ret] = snp_save_session(doch,event)
%% 
%  
%  file:   snp_save_session.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.06.29. Wednesday, 10:51:28
%

%% Requesting filename

persist = pcz_persist;

prompt = {['Save session as: ''' persist.date '-<fname>''']};
dlg_title = 'Save session';
num_lines = 1;
defaultans = {''};
answer = inputDialog(prompt,dlg_title,num_lines,defaultans);

if isempty(answer), return, end

answer = char(answer(1));

if isempty(answer)
    session_name = persist.date;
else
    session_name = [ persist.date '-' answer ];
end 

session_name = [ persist.session '/' session_name '.txt' ];

%% Save session

es = com.mathworks.mlservices.MLEditorServices;
editors = es.getEditorApplication.getOpenEditors;

fileID = fopen(session_name,'w');

for i = 0:editors.size()-1
    fname = editors.get(i).getLongName;
    
    fprintf(fileID, '%s\n', char(fname));
end

fclose(fileID);



end