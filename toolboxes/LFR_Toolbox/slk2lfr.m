%%
% SLK2LFR                From Simulink interconnections to LFR model
% ------------------------------------------------------------------
%
% This function  converts any Simulink  diagram  describing a linear 
% terconnection of several LFR objects  with LTI plants to  a single
% LFR object. 
% 
% The Simukink diagram is to be built from elements of the dedicated
% library LFRLIB.
%
% CALL
% sysLFR = slk2lfr(filename);
%
% INPUT ARGUMENT
% - filename: name of the Simulink file without mdl extension.
%        
% OUTPUT ARGUMENT
% - sysLFR: global LFR object.
