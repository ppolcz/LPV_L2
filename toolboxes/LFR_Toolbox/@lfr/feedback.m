% FEEDBACK     - Feedback for LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Computes a feedback loop involving one or two lfr-objects.
%
% SYNOPSIS
% sys3 = feedback(sys1,sys2[,sign[,nin[,nout]]]);
%
% DESCRIPTION
%
%      sys3 <=>      -->+---| sys1 |---|--->
%                       ^(-/+)         |
%                       |              |
%                       `---| sys2 |<--'
%
%
% INPUT ARGUMENTS
% sys1,sys2  One  is  an lfr-object, the  other  one can also be a
%            constant matrix or an LTI system.
% sign       -1 or +1.  By  default sign = -1  which means a minus
%            at the output of sys2.
% nin        specifies  the  inputs  (and  their ordering) of sys1
%            that are connected to the outputs of sys2
% nout       specifies  the  outputs  (and their ordering) of sys1
%            that are connected to the inputs of sys2
%
% OUTPUT ARGUMENT
% sys3      Resulting lfr-object
%
% See also lfr/minus, lfr/uminus, lfr/mrdivide, lfr/plus,
%          lfr/ctranspose, lfr/mtimes
%#----------------------------------------------------------------
