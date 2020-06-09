classdef he
%% He
%  
%  File: He.m
%  Directory: 7_ftools/utilities
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. September 01. (2019a)
%  Major review on 2020. May 22. (2019b)
%

%%


methods
    
    % 2020.05.22. (május 22, péntek), 00:25
    function [He,O,I] = he
        O = @(varargin) zeros(varargin{:});
        I = @(varargin) eye(varargin{:});
    end
    
    function varargout = subsref(~,S)
        varargout{1} = S(1).subs{1} + S(1).subs{1}';
    end
    
    % Ez mekkora hack!!
    function n = numel(varargin)
        n = 1;
    end

end


end
