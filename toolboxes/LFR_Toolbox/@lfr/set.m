% SET          - Changes parameter names, bounds, fields within LFR-object
%-----------------------------------------------------------------
% SYNOPSIS
% [sys = ]set(sys,par_name,new_par_name)
% [sys = ]set(sys,par_name,bound,btype)
% [sys = ]set(sys,fieldname,value)
%
% INPUT ARGUMENTS
% sys             lfr-object
% par_name        cell-array of parameter names to be changed
% new_par_name    cell-array of new parameter names
%
% bound           cell-array of new bound informations like:
%                                   - bound = {[min,max[,nom]]}
%                                   - bound = {[sec1,sec2]}
%                                   - bound = {ltisys(a,b,c,d)}
% btype           cell-array of bound-types like:
%                                   - btype = {'minmax'}
%                                   - btype = {'sec'}
%                                   - btype = {'freq'}
%                                   - btype = {'disc'}
%
% fieldname        name of field, e.g. 'a', 'b' ...
% value            value of sys.fieldname
%
% OUTPUT ARGUMENTS
% sys             lfr-object with changed field or parameter names
%
% See also lfr
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b c
%    sys = a+b+c;
%    size(sys);
%    par_name = {'a','c'};
%    new_par_name = {'anew','cnew'};
%    sys1 = set(sys,par_name,new_par_name);
%    size(sys1);
%    sys2 = set(sys1,new_par_name,{[2,7],[-6,4]},{'minmax','minmax'});
%    size(sys2);
%    sys3 = set(sys1,new_par_name,{[2,7,3],[4-6*i,4]},{'minmax','disc'});
%    size(sys3);
