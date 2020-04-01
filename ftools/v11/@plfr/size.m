function [q,m] = size(pLFR,dim)
%% size
%  
%  File: size.m
%  Directory: 7_ftools/ftools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 26. (2019b)
%

%%

    % disp 'plfr:size'
    % FIGYELEM: NAGYON SOK MINDEN MEGHIVJA EZT

    narginchk(1,2);

    [q,m] = size(pLFR.lfrtbx_obj);
    
    if nargin < 2
        nargoutchk(0,2);

        if nargout < 2
            q = [q m];
        end
    else
        nargoutchk(0,1);

        if dim == 2
            q = m;
        end
    end
end

