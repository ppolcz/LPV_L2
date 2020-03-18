classdef he
%% He
%  
%  File: He.m
%  Directory: 7_ftools/utilities
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. September 01. (2019a)
%

%%


methods
    
    function varargout = subsref(~,S)
        varargout{1} = S(1).subs{1} + S(1).subs{1}';
    end
    
    % Ez mekkora hack!!
    function n = numel(varargin)
        n = 1;
    end

end


end