% SIZE         - Size of an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Returns  the  number  of outputs, inputs  of the lfr-object sys.
% Displays  lfr-object  information  on  screen if invoked without
% output arguments.
%
% SYNOPSIS
% [no,ni,ns,nc,rd,cd] = size(sys)
%                     no : number of outputs
%                     ni : number of inputs
%                     ns : number of states
%                     nc : dimension of ConstBlock
%                     rd : row dimension of uncertainty matrix
%                          (without ConstBlock)
%                     cd : column dimension of uncertainty matrix
%                          (without ConstBlock)
%
% outputs = size(sys,1)        returns just the number of outputs
% inputs  = size(sys,2)        returns just the number of inputs
% states  = size(sys,'states') number of continuous/discrete states
% order   = size(sys,'order')  dimension of uncertainty block
%                              (without ConstBlock)
% npar    = size(sys,'p_name') dimension of uncertainty block 'p_name'
%
% size(sys)              displays information about inputs/outputs,
%                        states and uncertainty blocks in DELTA
%
% See also lfr, lfrdata, convert2lfr
