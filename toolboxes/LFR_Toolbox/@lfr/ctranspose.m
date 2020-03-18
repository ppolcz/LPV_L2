% CTRANSPOSE   - Conjugate transpose of an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Computes  the  conjugate  transpose of a lfr-object. Invoked by
% sys2 = sys1'
%
% REMARK
% For  dynamic  objects,  the conjugation is replaced by a transf-
% ormation of the form  (A,B,C,D) -> (-A',-C',B',D'). The function
% lfr/ctranspose.m  is  not  relevant  for  objects  with complex
% uncertainties.
%
% SYNOPSIS
% sys2 = sys1';
%
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs x y z
%    A = [x+2 1/y x;1 2 x;0 x 0];
%    B = [3*y 2; y 1;1/x 3];
%    C = [2 4 x];
%    D = [x y];
%    sys1 = abcd2lfr([A B;C D],3);
%    sys2 = sys1';
%
% % Verification
%    sys3 = abcd2lfr([-A' -C';B' D'],3);
%    distlfr(sys2,sys3)
