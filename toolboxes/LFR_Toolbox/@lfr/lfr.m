% LFR          - Core function for LFR-object generation
%-----------------------------------------------------------------
% PURPOSE
% Core function for lfr-object generation
%
% SYNOPSIS
% syslfr = lfr(m11,m12,m21,m22,blk)
% syslfr = lfr(Q)
% syslfr = lfr(a,b,c,d,'c')
% syslfr = lfr(a,b,c,d,'d')
% syslfr = lfr(NAME,TYPE,DIMS,BOUND,BTYPE)
%
% INPUTS
% m11,m12,m21,m22,blk: matrices of lfr-object
%
% Q: ss-object, pck-system or constant-matrix
%
% a,b,c,d: system matrices of a continuous('c')/discrete('d')-time
%          system
%
% Realization of lfr-object with one uncertainty-block, where:
% NAME:      string with name of uncertainty block
% TYPE:      composite string defining the uncertainty properties:
%            - nature: 'lti'   -> linear time-invariant (Default)
%                      'ltv'   -> linear time-varying
%                      'nl'    -> arbitrary nonlinear
%                      'nlm'   -> nonlinear memoryless
%            - structure: 's' -> scalar block d*eye(DIMS) (Default)
%                         'f' -> full block
%            - value:     'r' -> real-valued        (Default)
%                         'c' -> complex-valued
%            (e.g. 'ltifc' means lti-full-complex-block)
% DIMS:      dimension of uncertainty block:
%              DIMS = [1,2] means a 1x2 uncertainty block or
%              DIMS = 2 means a 2x2 uncertainty block
% BOUND:     quantitative information about the uncertainty:
%            - minmax-bounded uncertainty:
%                 set BOUND = [min,max] or
%                     BOUND = max (-> min=-max for norm-bound)
%            - frequency weighted bound:
%                 BOUND is SISO system W for frequency-
%                 weighted bounds:
%                       | Delta (jw) | <  | W(jw) |
%                 (e.g. BOUND=ltisys(-1,1,1,0))
%            - uncertainties in a disc
%                 set BOUND = [center,radius]
%                             (center is complex, radius is real)
%            - sector-bounded uncertainty:
%                 set  BOUND = [a,b]  for uncertainty valued in
%                 the sector {a,b}
% BTYPE:     - 'minmax' -> min/max-values bound (Default)
%            - 'freq'   -> frequency dependent bound
%            - 'sec'    -> sector bounded
%            - 'disc'   -> uncertainties in a disc
%
%
% DESCRIPTION
% An lfr-object is a generalized system in which Dynamics and
% uncertainty variations are modelled by an artificial feedback.
% So,  in  addition  to  the  matrices  (a,b,c,d),  an
% additional matrix 'blk' characterizes the blocks of the artificial
% feedback. The blocknames '1/s', '1/z' and '1' are
% reserved for the integrator block 'I/s' (continuous-time systems),
% the delay block 'I/z' (discrete-time systems) and constant block 'I'
% respectively, within the uncertainty matrix.
%
%
% 'blk' is a structure with
%
%  blk.desc : matrix  where   each  column  describes  one diagonal
%             block of the artificial feedback (uncertainty matrix).
%             Each column has 6 entries where
%
%             Entry no.      value    description
%
%                1            m       row dimension
%                2            n       column dimension
%                3           1/0      real(1)/complex(0)
%                4           1/0      scalar(1)/full(0)
%                5           1/0      linear(1)/nonlinear(0)
%                6           1/0      time-invariant(1)(memoryless
%                                     if nonlinear)/time-varying(0)
%              7:end                  bound information
%
%  blk.names: cell array including the names of the diagonal
%             blocks within the artificial feedback matrix
