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
%       2 1     % row-dimensions of blocks
%       2 1     % column-dimensions of blocks
%       1 1     % real(1) / complex(0) block types
%       1 1     % scalar(1) / full(0) block types
%       1 1     % linear(1) / nonlinear (0) block types
%       1 1     % time-inv.(1) / time-var.(0) block types
%       1 1     % min/max(1) / sector(2) / freq. dependent(>2) bounds, vagy 0
%       2 2     % min/max(2) / sector(1) / freq. dependent(>2) bounds, vagy 0
%       -1 -1   % minumum values of bounds
%       1 1     % maximum values of bounds
%       0 0 ]; 
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

%%

properties (Constant = true)
end

properties (GetAccess = private, SetAccess = private)
end

properties (GetAccess = public, SetAccess = private)
    lfrtbx_obj
    
    % Ez egy uj megoldas, gyors
    delta_fh

    M
    A,B,C,D,Delta,I
    np,nu,ny,m1
    
    % ebben azok a szimbolikus valtozok vannak, amibe kivulrol be akarok
    % helyettesiteni [AMIT MODOSITANI LEHET]
    subsvars % pl. [x1 p1 p2 x2 x3 dp1 dp2 .... barmi ]

    bounds
end

properties (GetAccess = public, SetAccess = public)
end

properties (Dependent)
    blk, names, 
    
    % ebben az 1 is benne lehet (lenyegeben Delta blokk neve)
    vars % pl. [1 x1 x2 p1 p2] 
    
    % ebben csak a szimbolikus valtozok lehetnek (amibe bele lehet
    % helyettesiteni
    symvars % pl. [x1 x2 p1 p2]
    
    % ebben azok a szimbolikus valtozok vannak, amibe kivulrol be akarok
    % helyettesiteni [IGY KEZELEM, HOGY ESETLEG NINCS BEALLITVA]
    % subsvars % pl. [x1 p1 p2 x2 x3 dp1 dp2 .... barmi ]

    desc
end

methods
    function blk = get.blk(pLFR)
        blk = pLFR.lfrtbx_obj.blk;
    end

    function names = get.names(pLFR)
        names = pLFR.lfrtbx_obj.blk.names;
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
            symvars = sym(pLFR.names(~strcmp('1',pLFR.names)));
        end
    end
    
    function desc = get.desc(pLFR)
        desc = pLFR.lfrtbx_obj.blk.desc;
    end
    
    % function subsvars = get.subsvars(pLFR)
    %     if isempty(pLFR.subsvars_)
    %         subsvars = pLFR.symvars;
    %     else
    %         subsvars = pLFR.subsvars_;
    %     end
    % end
end

methods (Static)
        
end

methods (Access = public)

    function pLFR = plfr(varargin) %(arg1,arg2,arg3,arg4,arg5,arg6)

        NrIn = nargin;
        
        if NrIn > 0 && isa(varargin{1},'plfr')
            pLFR = varargin{1};
            return
        end
        
        % Make everything EMPTY
        for prop = properties(pLFR).'
            try
                pLFR.(prop{1}) = '[EMPTY]';
            catch e
            end
        end
        
        % plfr(syslfr,subsvars)
        % 2019.11.20. (november 20, szerda), 15:30
        req_subsvars = [];
        if NrIn == 2 && isa(varargin{1},'lfr') ...
                && (isa(varargin{2},'sym') || isa(varargin{2},'lfr') || isempty(varargin{2}))
            if isa(varargin{2},'lfr')
                req_subsvars = sym(plfr(varargin{2}));
            else
                req_subsvars = varargin{2};
            end
            NrIn = 1;
        end            
        
        % Initialize: A,B,C,D, (lfrtbx_obj or (Delta and bounds))
        switch NrIn
            
            case 1 % (syslfr)
                pLFR.lfrtbx_obj = varargin{1};
                [pLFR.D,pLFR.C,pLFR.B,pLFR.A] = lfrdata(pLFR.lfrtbx_obj);

                
            case 2 % (M,blk)
                M = varargin{1};
                blk = varargin{2};
                m1 = sum(blk.desc(1,:));

                [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = pcz_split_matrix(M, [Inf m1], [Inf m1]);
                pLFR.lfrtbx_obj = lfr(pLFR.D,pLFR.C,pLFR.B,pLFR.A,blk);
                
                
            case 3 % (M,Delta,bounds)
                M = varargin{1};
                pLFR.Delta = varargin{2};
                pLFR.bounds = varargin{3};
                m1 = size(pLFR.Delta,1);

                [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = pcz_split_matrix(M, [Inf m1], [Inf m1]);
                
                
            case 5 % (A,B,C,D,blk)
                pLFR.lfrtbx_obj = lfr(varargin{[4,3,2,1,5]});
                [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = deal(varargin{1:4});
                
                
            case 6 % (A,B,C,D,5:Delta,6:bounds)
                [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = deal(varargin{1:4});
                pLFR.Delta = varargin{5};
                pLFR.bounds = varargin{6};
                
                if isempty(pLFR.bounds)
                    pLFR.bounds = ones(size(symvar(pLFR.Delta)))' * [-1 1];
                end

            % This case is more special compared to the others.
            case 7 % (A,B,C,D,5:Delta,6:[],7:blk)
                [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = deal(varargin{1:4});
                pLFR.Delta = varargin{5};
                pLFR.bounds = varargin{7}.desc(9:10,:)';
                pLFR = pLFR.generateLFR(varargin{7});
                
        end
        
        % Generate (Delta and bounds) or lfrtbx_obj
        switch NrIn
            
            case {1,2,5} % lfrtbx_obj --> (Delta, bounds)
                pLFR.bounds = pLFR.desc(9:10,:)';
                pLFR = pLFR.generateDelta;
                
            case {3,6} % (Delta, bounds) --> (lfrtbx_obj, bounds -- update)
                pLFR = pLFR.generateLFR;
                
        end
        
        % 2019.11.20. (november 20, szerda), 15:34
        if isempty(req_subsvars)
            req_subsvars = pLFR.symvars;
        end
        
        pLFR = pLFR...
            .set_vars(req_subsvars(:))... 2019.11.20. (november 20, szerda), 15:34
            .generateDimensions...
            .generateOtherFields...
            .test_properties;
    end

    % Overloading operators (), {} and .
    varargout = subsref(pLFR, S)
    
    varargout = val(plfr,varargin)

    pLFRT = transpose(pLFR)
    pLFRT = ctranspose(pLFR)
    ret = mtimes(X,Y);
    
    disp(pLFR);

    function H = plus(F,G)
        if isa(F,'plfr')
            F = F.lfrtbx_obj;
        end

        if isa(G,'plfr')
            G = G.lfrtbx_obj;
        end
        H = plfr(F + G);
    end
    
    function H = minus(F,G)
        if isa(F,'plfr')
            F = F.lfrtbx_obj;
        end

        if isa(G,'plfr')
            G = G.lfrtbx_obj;
        end
        H = plfr(F - G);
    end
    
    % 2019.12.09. (december  9, hétfő), 01:07
    function H = mldivide(F,G)
        if isa(F,'plfr')
            F = F.lfrtbx_obj;
        end

        if isa(G,'plfr')
            G = G.lfrtbx_obj;
        end
        H = plfr(mldivide(F,G));
    end
    
    % 2019.12.09. (december  9, hétfő), 01:07
    function H = mrdivide(F,G)
        if isa(F,'plfr')
            F = F.lfrtbx_obj;
        end

        if isa(G,'plfr')
            G = G.lfrtbx_obj;
        end
        H = plfr(mrdivide(F,G));
    end
    
    % 2019.12.09. (december  9, hétfő), 00:56
    function syslfr = lfr(pLFR)
        syslfr = pLFR.lfrtbx_obj;
    end
    
    % 2019.11.21. (november 21, csütörtök), 10:29
    function H = vertcat(varargin)
        H = plfr.cat__(@vertcat,varargin{:});
    end

    % 2019.11.21. (november 21, csütörtök), 10:29
    function H = horzcat(varargin)
        H = plfr.cat__(@horzcat,varargin{:});
    end

    function pLFR = set_vars(pLFR,vars_)
        
        nem_szerepel = setdiff(pLFR.symvars,vars_);

        if ~isempty(nem_szerepel)
            error('plfr: Variables `%s` must be included!', ...
                strjoin(cellfun(@(x) {char(x)}, num2cell(nem_szerepel)),', '))
        end
        
        pLFR.subsvars = vars_(:);
                
        delta = diag(pLFR.Delta);
        
        if isempty(delta)
            pLFR.delta_fh = @(xp) zeros(0,size(xp,2));
        else
            ZERO = sym('ZERO');
            delta(delta == 1) = 1 + ZERO;

            pLFR.delta_fh = matlabFunction(delta,'vars',[ pLFR.subsvars(:) ; ZERO ]);
        end
    end
    
    function pLFR = generateDimensions(pLFR)

        % LFR dimensions (output, input)
        [pLFR.ny,pLFR.nu] = size(pLFR.A);

        % LFR dimensions (Delta block)
        pLFR.m1 = size(pLFR.Delta,1);

        % LFR dimensions (nr. of parameters)
        pLFR.np = numel(pLFR.symvars);        
    end
    
    % TODO: nem lenne-e jobb dependent var-ket kezelni?
    function pLFR = generateOtherFields(pLFR)
        pLFR.M = [pLFR.A,pLFR.B;pLFR.C,pLFR.D];
        pLFR.I = eye(pLFR.m1);        
    end
    
    % Requires: valid lfrtbx_obj
    function pLFR = generateDelta(pLFR)
        s = numel(pLFR.names);
        c = cell(s, 1);
        for k = 1:s
            c{k} = sym(pLFR.names{k}) * eye(pLFR.desc(1:2,k)');
        end
        pLFR.Delta = blkdiag(c{:});
    end

    function pLFR = generateLFR(pLFR,blk)

        % Elotte_Gx = pLFR.A + pLFR.B/(I - pLFR.Delta*pLFR.D)*pLFR.Delta*pLFR.C;

        [delta,sigma] = sort(diag(pLFR.Delta));

        pLFR.Delta = diag(delta);

        Is = pcz_permat(sigma);

        pLFR.A = pLFR.A;
        pLFR.B = pLFR.B*Is;
        pLFR.C = Is'*pLFR.C;
        pLFR.D = Is'*pLFR.D*Is;
        pLFR.Delta = Is' * pLFR.Delta * Is;

        % Utana_Gx = A + B/(I - Delta*D)*Delta*C;
        % pcz_symzero_report(tGx - Gx, 'new G(x) = G(x) before permutation');

        [p,cumnr] = unique(delta);

        % If the bounds does not contain the first 0 0 values for the
        % constant part of the Delta block
        if p(1) == 1 && size(pLFR.bounds,1) == numel(p)-1
            pLFR.bounds = [
                0 0
                pLFR.bounds
                ];
        end
        
        dim = diff([ cumnr ; size(pLFR.Delta,1)+1 ])';
        
        % dim
        % pLFR.bounds'

        lfrtbx_blk.names = cellfun(@(pi) {char(pi)}, num2cell(p.'));

        if nargin > 1
        
            % We have already a desc matrix, but we need to update it.
            
            % lfrtbx_blk.names
            % blk.names
            
            idx = zeros(1,numel(lfrtbx_blk.names));
            for i = 1:numel(lfrtbx_blk.names)
                idx(i) = find(strcmp(blk.names,lfrtbx_blk.names{i}));
            end
            
            lfrtbx_blk.desc = blk.desc(:,idx);
            
            lfrtbx_blk.desc(1,:) = dim;
            lfrtbx_blk.desc(2,:) = dim;
                        
        else

            % Delta interface for LFR Toolbox
            lfrtbx_blk.desc = [
                % row dimension
                dim
                % column dimension
                dim
                % real(1)/complex(0)
                dim*0 + 1
                % scalar(1)/full(0)
                dim*0 + 1
                % linear(1)/nonlinear(0)
                dim*0 + 1
                % time-invariant(1)(memoryless if nonlinear)/time-varying(0)
                dim*0 + 1
                % bound informations (given as intervals: 1,2)
                dim*0 + 1
                dim*0 + 2
                % b.i.: intervals
                pLFR.bounds'
                % b.i.: nominal values
                [0.5 0.5]*pLFR.bounds'
                ];

        end
            
        pLFR.lfrtbx_obj = lfr(pLFR.D,pLFR.C,pLFR.B,pLFR.A,lfrtbx_blk);
    end

    function pLFR = test_properties(pLFR)

        failed = false;
        
        for prop = properties(pLFR).'
            if strcmp(pLFR.(prop{1}),'[EMPTY]')
                fprintf('Property %s is empty\n', prop{1});
                failed = true;
            end
        end
        
        if failed
            pcz_dispFunctionStackTrace('first',0,'last',0);
        end

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
    
    % 2019.11.20. (november 20, szerda), 14:46
    function [pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D)
        PI_sym = [
            eye(pLFR.nu)
            (pLFR.I-pLFR.Delta*D)\pLFR.Delta*C 
            ];
        
        if ~isnumeric(PI_sym)
            PI_sym = simplify(PI_sym);
        end
        
        pLFR_sym = [A B]*PI_sym;
        
        if ~isnumeric(pLFR_sym)
            pLFR_sym = simplify(pLFR_sym);
        end        
    end
    
    % Atirva: 2019.11.20. (november 20, szerda), 14:46
    function [pLFR_sym, PI_sym] = sym(pLFR)
        [A,B,C,D] = data(pLFR); %#ok<PROP>        
        [pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D); %#ok<PROP>
    end
    
    % Created: 2019.11.20. (november 20, szerda), 14:46
    function [pLFR_sym, PI_sym, fancy, s] = symr(pLFR,~)
        [A,B,C,D] = data(pLFR); %#ok<PROPLC>
        
        r = @(M) pcz_find_recdec(M,'maxden',1000,'maxit',1,'tol',1e-10, 'structout', true);

        [M,s] = r(pLFR.M);  %#ok<PROPLC>
        
        fancy = s.good;
        
        if nargin > 1 && ~fancy
            pLFR_sym = [];
            PI_sym = [];
            return
        end
            
        if fancy
            [A,B,C,D] = pcz_split_matrix(M,[Inf pLFR.m1],[Inf pLFR.m1]);  %#ok<PROPLC>
        end

        [pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D);  %#ok<PROPLC>
    end
    
    function PI = generatePI(pLFR)
        s1 = pLFR.ny;
        s2 = pLFR.nu;
        m = s2+pLFR.m1;

        [M11,M12] = pcz_split_matrix(eye(m),m,[s2 pLFR.m1]);

        PI = plfr(M11,M12,pLFR.C,pLFR.D,pLFR.blk);
        PI = PI.set_vars(pLFR.subsvars(:));
    end
    
    function PI1 = generatePI1(pLFR)
        s1 = pLFR.ny;
        s2 = pLFR.nu;
        m = s2+pLFR.m1;

        [~,~,M11,M12] = pcz_split_matrix(eye(m),[s2 pLFR.m1],[s2 pLFR.m1]);

        PI1 = plfr(M11,M12,pLFR.C,pLFR.D,pLFR.blk);
        PI1 = PI1.set_vars(pLFR.subsvars(:));
    end
    
    function pLFR = minlfr(pLFR)
        pLFR = plfr(minlfr(pLFR.lfrtbx_obj),pLFR.subsvars(:));
    end
    
    function dLFR = diff(pLFR,p_cell,dp_cell,subsvars)
        
        dLFR = minlfr(pLFR*0);
        
        for i = 1:numel(p_cell)
            try
                dLFR = dLFR + diff(pLFR.lfrtbx_obj,p_cell{i}.blk.names{1}) * dp_cell{i};
            catch e
                if strfind(e.message, 'does not exist in lfr-object')
                    continue
                else
                    getReport(e)
                    break;
                end
            end
        end
        
        dp_plfr = plfr(vertcat(dp_cell{:}));
        
        if nargin < 4
            subsvars = plfr.var_helper__(pLFR,dp_plfr);
        end
        
        dLFR = dLFR.set_vars(subsvars(:));                
    end
    
    N = paffmat(pLFR,xi);
        
end

methods (Static, Access = private)
%% private    

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

% 2019.11.21. (november 21, csütörtök), 10:29
function H = cat__(catfn, varargin)
    if nargin < 2
        H = [];
    else
        varargin_lfr = cellfun(@(l) {plfr.to_lfr(l)}, varargin);
        varargin_plfr = cellfun(@(l) {plfr.to_plfr(lfr(l))}, varargin);
        H = plfr(catfn(varargin_lfr{:}));
        H = H.set_vars(plfr.var_helper__(varargin_plfr{:}));
    end
end


end

end


