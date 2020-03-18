function pLFRT = transpose(pLFR)
%
%  File: transpose.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2019. March 29.

pLFRT = plfr(pLFR.lfrtbx_obj.',pLFR.subsvars);

end

