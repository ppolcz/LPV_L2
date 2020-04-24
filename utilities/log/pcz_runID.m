function [ret, A_runID_dates] = pcz_runID(msg1, msg2, varargin)
%%
%  file:   pcz_runID.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.05.11. Wednesday, 15:08:37
%  Major review on 2020. April 22. (2019b)
%

A_MAT_PATH = getenv('A_MAT_PATH');

try  
    load(A_MAT_PATH);    
catch
    %%
    [path,bname,ext] = fileparts(A_MAT_PATH);
    [status,msg,msgID] = mkdir(path);
    fprintf('I attempted to create folder %s.\n', path);
    fprintf('Status: %d\n', status);
    fprintf('Message: %s\n', msg);
    fprintf('Message ID: %s\n', msgID);
    A_runID_dates = cell2table(cell(0,4),'VariableNames',{'RunID', 'Date', 'Comment', 'Reserved'});
    A_runID = -1;
    
    fprintf('I will create registry `%s%s'' for logging reasons.\n', bname, ext);
end

if nargin == 1 && strcmp(msg1,'-l')

    ret = A_runID_dates;
    
else
    
    if nargin <= 1
        msg2 = '';
    end

    if nargin <= 0
        msg1 = '';
    end

    % Increment A_runID
    A_runID = A_runID + 1;
    ret = sprintf('%04d',A_runID);

    % Register new runID
    A_runID_dates = [ A_runID_dates ; {ret , pcz_fancyDate , msg1 , msg2} ];

    % Save variables
    save(A_MAT_PATH, 'A_runID', 'A_runID_dates')

end

end