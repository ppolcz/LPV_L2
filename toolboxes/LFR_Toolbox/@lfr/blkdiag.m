% BLKDIAG      - Diagonal concatenation of LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Diagonal concatenation of lfr-objects. This function can be used
% with lfr-objects having no input/output.
%                      | sys1   0  |
%                      |  0   sys2 |
%
% SYNOPSIS
% sys = blkdiag(sys1,sys2);
% sys = blkdiag(sys1,sys2,sys3,...);
%
% See also lfr/horzcat lfr/vertcat
%#----------------------------------------------------------------
% % EXAMPLE
%    sys1 = rlfr(3,2,3,1,0,1,'x');
%    sys2 = rlfr(2,1,3,0,1,2,'y');
%    distlfr([sys1 zeros(2,3);zeros(1,3) sys2],blkdiag(sys1,sys2))
