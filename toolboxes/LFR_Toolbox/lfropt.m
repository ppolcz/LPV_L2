% LFROPT       - Function for systematic order reduction
%-----------------------------------------------------------------
% PURPOSE
% This  function  defines global parameters that are used by other
% functions  of  the toolbox for applying systematically minlfr or
% minlfr1 with a given tolerance argument.
%
% SYNOPSIS
% lfropt             <- initializes global parameters if necessary
% lfropt('nD'[,tol]) <- applies minlfr after all operations
% lfropt('1D'[,tol]) <- applies minlfr1 after all operations
% lfropt('info')     <- for information
% lfropt('default')  <- default (no systematic order reduction)
%#----------------------------------------------------------------
% % EXAMPLES
%    clear global
%    lfropt('i')
%    lfropt
%    lfropt('i')
%    lfropt('1',0.0001)
%    lfropt('i')
%    lfropt('n',0.0001)
%    lfropt('i')
