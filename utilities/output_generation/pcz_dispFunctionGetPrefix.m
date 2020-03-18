function [prefix] = pcz_dispFunctionGetPrefix
%% pcz_dispFunctionGetPrefix
%  
%  File: pcz_dispFunctionGetPrefix.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. April 30.
%

%%

global SCOPE_DEPTH

% [ST,I] = dbstack;
% 
% for i = 2:SCOPE_DEPTH
%     fprintf('│   ')
% end

% if numel(ST) > I    
%     if ~isempty(msg)
%         disp(['│   - ' msg])
%     else
%         disp '│   '
%     end
% else
%     disp(['- ' msg ])
% end

prefix = '';
if SCOPE_DEPTH >= 1
    tab = '│   ';
    prefix = repmat(tab,[1 SCOPE_DEPTH]);
end


end