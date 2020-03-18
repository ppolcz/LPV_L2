% REORDERLFR   - Re-ordering lfr uncertainty blocks
%-----------------------------------------------------------------
% PURPOSE
% This function is always called at the end of an LFT-realization
% and performs lexicographical ordering of the uncertainty blocks
% in the Delta-Matrix. Furthermore, repeated  scalar  blocks  are
% merged  to  one  block  in  the  uncertainty  description.  The
% special block with name "ConstBlock" is put  on  first position
% of the Delta-Matrix. The blocks "Int" (or "Delay" for  discrete
% time systems) is put on the second position of the Delta-Matrix
% (or on the first position, if no "ConstBlock" exists).
% Then Nonlinear  and Time-varying blocks  are put in  third  and 
% fourth poistions (or second and third or first and second).
%
% SYNOPSIS
% [a,b,c,e,blk] = reorderlfr(a1,b1,c1,e1,blk1)
%
% INPUT ARGUMENTS
% a1,b1,c1,e1,blk1     lfr-object matrices
%                      (Note: d-matrix neglected as it is not
%                       changed by reorderlfr!)
%
% OUTPUT ARGUMENT
% a,b,c,e,blk          ordered lfr-object matrices
%#----------------------------------------------------------------
