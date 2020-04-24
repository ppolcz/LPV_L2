classdef blk
%% blk
%  
%  File: blk.m
%  Directory: 7_ftools/ftools/v12/@blk
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. April 10. (2019b)
%

%%
properties (Constant = true)
end

properties (GetAccess = private, SetAccess = private)
    
    % Array of struct objects. E.g.:
    %              Name: 'x1'
    %              Type: 'PAR'
    %              Size: [13 13]
    %          NomValue: 0
    %            Bounds: [-1 1]
    %        RateBounds: [-1 1]
    %          Measured: 'no'
    %     Normalization: []
    %              Misc: []    
    s 
end

methods (Access = public)

    
    varargout = subsref(D, S)
end

methods (Access = private)
end

end
