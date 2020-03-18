% LFRT         - List of functions of the LFRT on 10-10-2005
%
% - Overloaded functions
%     display        - show contents of lfr-object.
%     size           - show size information of lfr-object.
%     isempty        - check if lfr-object is empty.
%     get            - get fields of an lfr-object.
%     set            - change block-names, block-bounds or object-fields.
%     plus           - addition of lfr-objects.
%     minus          - subtraction of lfr-objects.
%     uminus         - sign change of lfr-object.
%     mtimes         - multiplication of lfr-objects.
%     mpower         - repeated multiplication of lfr-object.
%     mrdivide       - division of lfr-objects.
%     inv            - inversion of lfr-object.
%     horzcat        - horizontal concatenation of lfr-objects.
%     vertcat        - vertical concatenation of lfr-objects.
%     append         - block diagonal concatenation
%     blkdiag        - block diagonal concatenation
%     transp         - transposition of lfr-objects.
%     subsref        - extract matrices of lfr-object.
%     subsasgn       - assignment function for lfr-objects.
%     feedback       - feedback interconnection of lfr-objects.
%     dcgain         - steady state gain
%     eval           - evaluation from values in the workspace
%     diff           - differentiation
%     eig            - nominal poles of a dynamic lfr-object
%
% - lfr-object generation
%     lfr            - core function of the toolbox
%     lfrs           - realize elementary real or complex lfr-objects.
%     rlfr           - generates random lfrs
%     bnds2lfr       - from min/max bounds to lfr-object
%     sym2lfr        - from symbolic to lfr-object
%     symtreed       - elementary structured tree decomposition
%     aeblkchk       - check well-definedness of lfr-object.
%
%  - Conversions
%     lfr            - convert various objects to lfr-objects
%     abcd2lfr       - from state-space to input-output form
%     lfr2abcd       - converse
%     lf2lfr         - lfr-object from left fractional factorization.
%     rf2lfr         - lfr-object from right fractional factorization.
%     lfr2rob        - transforms an lfr-object to umat or uss-object
%
%  - lfr realization
%     sym2lfr        - symbolic expression converted to lfr
%     symtreed       - elementary structured tree decomposition
%     lfrdata        - retrieves realization matrices from an lfr-object
%     gmorton        - generalized Morton realization
%
%  - Gridding and interpolation
%     data2lfr       - interpolation on an lfr basis
%     plotlfr        - plots entries of an lfr object on a gridding
%
%  - Order reduction
%     minlfr         - reduction based on n-D Kalman decomposition
%     minlfr1        - sequencial 1-D odrer reduction
%     reduclfr       - approximation based on the above tools
%     symtreed       - elementary structured tree decomposition
%
%  - Manipulation of the delta-block
%     eval           - evaluates an lfr from values in the workspace
%     uplft          - closes (partially) the delta-loop
%     normalizelfr   - normalizes parametric blocks in lfr-object
%     unnormalize    - converse
%     actualval      - retrieves actual values from normalized ones
%     starplfr       - star product
%
%  - Dynamic systems
%     feedback       - closes the feedback loop
%     dcgain         - dc-gain
%     eig            - nominal poles
%
%  - Mu analysis
%     lfr2mustab     - from lfr-object to input arguments of mustab
%     lfr2mubnd      - from lfr-object to input arguments of mubnd
%     lfr2mu         - from lfr-object to input arguments of mu
%     lfr2mussv      - from lfr-object to input arguments of mussv
%     ns_rad         - non-singularity radius
%     wp_rad         - well-posedness radius
%     min_max        - min and max values of a 1-by-1 real lfr-object
%
%  - Distances
%     distlfr        - calculates a lower bound of the distance between two lfrs
%     udistlfr       - upper bound of the distance (some restrictions apply)
