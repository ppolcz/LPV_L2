function [q,m] = size(N,dim)
%%
%  File: size.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019.09.08. (szeptember  8, vasárnap), 13:07 [ATIRVA - OK]

%%
    % disp 'PGenAffineMatrix:size'
    % FIGYELEM: NAGYON SOK MINDEN MEGHIVJA EZT

    narginchk(1,2);

    [q,msp1] = size(N.Theta);
    m = msp1 / numel(N.channels);

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

