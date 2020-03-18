% LFRS         - Generation of real/complex elementary LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Shortcut  for  defining  real or complex elementary lfr-objects.
% The  names  'Int'  and  'Delay' are reserved for the Integration
% (1/s) and Delay (1/z) operator respectively.
%
% SYNOPSIS
% lfrs a b [real] [min max [nom]]  generates real elementary
%                                  lfr objects a,b
%                                  (the keyword real is optional)
% lfrs a b complex [center radius] generates complex elementary
%                                  lfr objects a,b
%
% INPUT ARGUMENTS
% a b         Strings:  names  of  the  elementary real or complex
%             lfr-objects   to  be  created. The  key word complex
%             determines if the parameters are complex or real.
% min         vector of minimum values (default: -1) for variables
%             (for  'Int'/'Delay' one may choose arbitrary values)
% max         vector  of maximum values (default: 1) for variables
%             (for  'Int'/'Delay' one may choose arbitrary values)
% nom         vector of nominal values (default: 0)
% center      vector of centers, entries can be complex (default: 0)
% radius      vector of radius, real entries (default: 1)
%
% NOTE
% If  the  vectors  min,  max,  nom, center or radius occur, their
% length  must  be equal to the amount of variables. Even if there
% is  a  single variable, the above vectors must be strings of the
% form [...] (don't forget the square brackets).
%
% OUTPUT ARGUMENTS
% Generated    elementary lfr-objects with name of input arguments.
%
% See also lfr
%-----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b c
%    lfrs d e f complex
%    lfrs a b c [1,2,3] [4,5,6]
%    lfrs a b c real [1,2,3] [4,5,6]
%    lfrs a b complex [-1,-2] [1,2]
%    lfrs a Int [50,0] [100,0]
