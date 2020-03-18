% STARPLFR     - Star product of LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Star product of LFR-objects
%
% SYNOPSIS
% sysout = starplfr(syslow,sysup);
%
% WARNING
% The  number of outputs and inputs of the upper LFR must be equal
% to  the  size the global delta-block of the lower LFR. Note that
% if  you  don't  want to replace ALL the Delta block of the lower
% LFR by the upper LFR, use the function 'eval'.
%
% INPUT ARGUMENTS
% syslow   Lower LFR
% sysup    Upper LFR
%
% OUTPUT ARGUMENT
% sysout   Returned LFR.
%
% See also lfr/eval, uplft
%#----------------------------------------------------------------
% % EXAMPLE
%     blklow = struct('names',{{'A'}},'desc',[6;6;0;0;1;1;1;2;-1;1]);
%     syslow = rlfr(2,3,blklow);   % the Delta block of sysout is 6-by-6
%     sysup  = rlfr(0,6,6,0,3);    % sysup is 6-by-6
%     sysout1 = starplfr(syslow,sysup);
%
% % Equivalent
%     A = sysup;
%     sysout2 = eval(syslow);
%     distlfr(sysout1,sysout2)
