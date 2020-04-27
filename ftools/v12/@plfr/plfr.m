classdef plfr
% plfr
%
%  File: plfr.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. February 19.
%
%
% Structure of the LFR Toolbox (source: https://www.researchgate.net/profile/Simon_Hecker/publication/4124334_Enhanced_LFR-toolbox_for_MATLAB/links/0046352c19a1b2f676000000/Enhanced-LFR-toolbox-for-MATLAB.pdf)
%
%   blk.names = {'p1','p2'};
%   blk.desc = [
%     1: 2  2  1     % row-dimensions of blocks
%     2: 2  2  1     % column-dimensions of blocks
%     3: 1  1  1     % real(1) / complex(0) block types
%     4: 1  1  1     % scalar(1) / full(0) block types
%     5: 1  1  1     % linear(1) / nonlinear (0) block types
%     6: 1  1  1     % time-inv.(1) / time-var.(0) block types
%     7: 0  1  1     % min/max(1) / sector(2) / freq. dependent(>2) bounds, vagy 0
%     8: 0  2  2     % min/max(2) / sector(1) / freq. dependent(>2) bounds, vagy 0
%     9: 0 -1 -1     % minumum values of bounds
%    10: 0  1  1     % maximum values of bounds
%    11: 0  0  0     % nominal values
%    -------------------------------------
%    12: 0 -1 -2     % minimum rate bound
%    13: 0  1  2     % maximum rate bound
%       ];
%
%
%
% Lower LFR: G(p) = A + B*(I - Delta*D)^{-1}*Delta*C
%
%  M = [A B ; C D]: (nu+m1)x(ny+m1)
%
%  A: nu x ny
%  B: nu x m1
%  C: m1 x ny
%  D: m1 x m1
%  I: m1 x m1
%
%  np nr. of parameter values (p \in \mathbb{R}^{np})
%
%    Assume that p1 \in [-1 1]
%                p2 \in [2 5]
%                p3 \in [-1 0]
%    then, bounds = [
%              -1 1
%               2 5
%              -1 0
%               ];
%
%    Note that if Delta containts a constant block (with ones in the
%    diagonal), then, bounds = [
%                         0 0
%                        -1 1
%                         2 5
%                        -1 0
%                         ];
%
% Preferred order of variables:
%
%  x1,x2,p1,p2,dp1,pd2
%

%%

properties (Constant = true)

    Preferred_Variable_Order = struct(...
        'x', 1, ...
        'p', 10, ...
        'dp', 20, ...
        'other', Inf ...
        );
    Parse_VarName = @(varname) regexp(varname,'(?<name>[a-zA-Z]+)(?<nr>[0-9_]*)','once','names');
end

properties (GetAccess = private, SetAccess = private)
end

properties (GetAccess = public, SetAccess = private)
    lfrtbx_obj % --> blk, desc, names, vars, symvars

    % Fast evaluation in multiple (numerical) points
    delta_fh

    % Flexible evaluation in a single symbolical point of different types
    % deltai_fh_cell

    M
    A,B,C,D,Delta,I
    nu,ny,m1

    % ebben azok a szimbolikus valtozok vannak, amibe kivulrol be akarok
    % helyettesiteni
    % pl. [x1 p1 p2 x2 x3 dp1 dp2 .... barmi ]
    subsvars % --> ns
end

properties (GetAccess = public, SetAccess = public)
end

properties (Dependent)
    blk, desc, names, bounds, rbounds

    varnames % pl. { 'x1' 'x2' 'p1' 'p2' }

    % ebben az 1 is benne lehet (lenyegeben Delta blokk neve)
    vars % pl. [1 x1 x2 p1 p2]

    % ebben csak a szimbolikus valtozok lehetnek (amibe bele lehet
    % helyettesiteni
    symvars % pl. [x1 x2 p1 p2]

    % Nr. of variables:
    % E.g. A(p1,p2,p3) = [ 1+p1+p2 , p2^2 ]
    np % number of variables (expected from outside): 3
    ns % actual number of variables                 : 2
end

methods
    function np = get.np(pLFR)
        np = numel(pLFR.subsvars);
    end

    function ns = get.ns(pLFR)
        ns = numel(pLFR.symvars);
    end

    function blk = get.blk(pLFR)
        blk = pLFR.lfrtbx_obj.blk;
    end

    function names = get.names(pLFR)
        names = pLFR.lfrtbx_obj.blk.names;
    end

    function varnames = get.varnames(pLFR)
        varnames = pLFR.names(~strcmp('1',pLFR.names));
    end

    function vars = get.vars(pLFR)
        if isempty(pLFR.names)
            vars = [];
        else
            vars = sym(pLFR.names).';
        end
    end

    function symvars = get.symvars(pLFR)
        if isempty(pLFR.names)
            symvars = [];
        else
            symvars = sym(pLFR.varnames);
        end
    end

    function desc = get.desc(pLFR)
        desc = pLFR.lfrtbx_obj.blk.desc;
    end

    function bounds = get.bounds(pLFR)
        bounds = pLFR.desc(9:10,sum(abs(pLFR.desc(7:11,:))) > 0)';
    end

    function rbounds = get.rbounds(pLFR)
        rbounds = pLFR.desc(12:13,sum(abs(pLFR.desc(7:11,:))) > 0)';
    end

end

methods (Static)
end

methods (Access = public)

    % Constructor
    function pLFR = plfr(varargin)
        pLFR = pLFR.init(varargin{:});
    end

end

methods (Access = private)

    % Constructor functions
    pLFR = init(pLFR,varargin)
    pLFR = generate_dimensions(pLFR)
	pLFR = generate_other_fields(pLFR)
	pLFR = generate_delta(pLFR)
	pLFR = generate_LFR(pLFR,blk)
    pLFR = test_properties(pLFR)

end

methods (Access = public)

    % Only for single value substitution: pLFR(1,2,3,4)
    varargout = subsref(pLFR, S)

    % Fast substitution for vector values: pLFR.val(grid_p1,grid_p2,...)
    varargout = val(plfr,varargin)

    [q,m] = size(pLFR,dim)
    disp(pLFR);

    % Modify substitution variables (external interface)
    pLFR = sort_vars(pLFR)
    pLFR = set_vars(pLFR,vars_)

    [pLFR_sym, PI_sym] = sym(pLFR)
    [pLFR_sym, PI_sym, fancy, s] = symr(pLFR,varargin)

    PI = generatePI(pLFR)
    PI1 = generatePI1(pLFR)

    dLFR = diff(pLFR,p_cell,dp_cell,subsvars)

    % Cast to PAffineMatrix object
    N = paffmat(pLFR,xi);

    % Cast to lfr object
    function syslfr = lfr(pLFR)
        syslfr = pLFR.lfrtbx_obj;
    end

    function [A,B,C,D,Delta,bounds,blk] = data(pLFR)
        A = pLFR.A;
        B = pLFR.B;
        C = pLFR.C;
        D = pLFR.D;
        Delta = pLFR.Delta;
        bounds = pLFR.bounds;
        blk = pLFR.blk;
    end

    % 2020.04.27. (április 27, hétfő), 12:26
    function data = export(pLFR)
        data.type = 'Lower LFR: M11 + M12 (I - Delta M22) Delta M22';
        data.M11 = pLFR.A;
        data.M12 = pLFR.B;
        data.M21 = pLFR.C;
        data.M22 = pLFR.D;
        data.Delta_str = sprintf('diag([%s])',strjoin(cellfun(@(s) {char(s)}, num2cell(diag(pLFR.Delta))),','));
        data.Delta_blk = cell2table([ num2cell(pLFR.blk.desc) {
            "Row  1. Row-dimensions of blocks"
            "Row  2. Column-dimensions of blocks"
            "Row  3. Real(1) / complex(0) block types"
            "Row  4. Scalar(1) / full(0) block types"
            "Row  5. Linear(1) / nonlinear (0) block types"
            "Row  6. Time-inv.(1) / time-var.(0) block types"
            "Row  7. Min/max(1) / sector(2) / freq. dependent(>2) bounds, 0 for constant block"
            "Row  8. Min/max(2) / sector(1) / freq. dependent(>2) bounds, 0 for constant block"
            "Row  9. Minumum values of bounds, 0 for constant block"
            "Row 10. Maximum values of bounds, 0 for constant block"
            "Row 11. Nominal values, 0 for constant block"
            "Row 12. Minimum rate bound"
            "Row 13. Maximum rate bound"
            } ], 'VariableNames', [ pLFR.blk.names 'LFR Toolbox block description']);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% Simplified operator wrappers
  %  2020.03.29. (március 29, vasárnap), 12:55

    function H = plus(varargin)
        H = plfr.operator__(@plus,varargin{:});
    end

    function H = minus(varargin)
        H = plfr.operator__(@minus,varargin{:});
    end

    function H = mtimes(varargin)
        H = plfr.operator__(@mtimes,varargin{:});
    end

    function H = mldivide(varargin)
        H = plfr.operator__(@mldivide,varargin{:});
    end

    function H = mrdivide(varargin)
        H = plfr.operator__(@mrdivide,varargin{:});
    end

    function H = vertcat(varargin)
        H = plfr.operator__(@vertcat,varargin{:});
    end

    function H = horzcat(varargin)
        H = plfr.operator__(@horzcat,varargin{:});
    end

    function H = ctranspose(G)
        H = plfr(G.lfrtbx_obj',G.subsvars(:));
    end

    function H = transpose(G)
        H = plfr(G.lfrtbx_obj.',G.subsvars(:));
    end

    function H = minlfr(G)
        H = plfr(minlfr(G.lfrtbx_obj),G.subsvars(:));
    end

  %
 %%% Simplified operators [END]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

methods (Access = private)

    [pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D)
    order_nr = ordernr__(pLFR,varname)

end


methods (Static, Access = private)

    function subsvars = var_helper__3(allvars)
        [~,ind] = unique(allvars);
        subsvars = allvars(sort(ind));
    end

    function subsvars = var_helper__2(varargin)
        varargin = cellfun(@(v){v(:)}, varargin);
        subsvars = plfr.var_helper__3(vertcat(varargin{:}));
    end

    function subsvars = var_helper__(varargin)
        allvars = cellfun(@(X) {X.subsvars(:)}, varargin);
        subsvars = plfr.var_helper__2(allvars{:});
    end

    function F = to_lfr(F)
        if isa(F,'plfr')
            F = F.lfrtbx_obj;
        end
    end

    function F = to_plfr(F,subsvars)
        if nargin < 2
            subsvars = [];
        end

        F = plfr(F, subsvars(:));
    end

    % Created:  2019.11.21. (november 21, csütörtök), 10:29
    % Modified: 2020.03.29. (március 29, vasárnap), 12:49
    function H = operator__(catfn, varargin)
        if nargin < 2
            H = [];
        else
            varargin_lfr = cellfun(@(l) {plfr.to_lfr(l)}, varargin);
            H = sort_vars(plfr(catfn(varargin_lfr{:})));
        end
    end


end

end


