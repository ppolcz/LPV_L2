function varargout = val(pLFR,varargin)
%%
%  File: val.m
%  Directory: 7_ftools/ftools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. September 05. (2019a)
%
%  1. case: varargin{1} = [ p11   p21   ...   pN1
%                           ...   ...   ...   ...
%                           p1np  p2np  ...   pNnp ]
%           assumption: pLFR is vertical!
%           output:     horzcat({:})
%
%  2. case: varargin{1} = [ p11  ...  p1np
%                           p21  ...  p2np
%                           ...  ...  ...
%                           pN1  ...  pNnp ]
%            assumption: pLFR is horizontal
%            output:     vertcat({:})
%
%  3. case: varargin{:} = { grid_p1 , ... , grid_pN }
%           assumption: pLFR is scalar!
%           output:     grid_pLFR

%%

CASE = 0;

% CASE = 1 or 2
if numel(varargin) == 1 && isnumeric(varargin{1})

    samples = varargin{1};

    CASE = find(size(samples) == pLFR.np);

    if isempty(CASE)
        error('Dimension mismach: np = %d, size(arg1) = [%d,%d].', pLFR.np, size(samples,1), size(samples,2))
    end

    if CASE == 2
        samples = samples';
    end

elseif numel(varargin) == pLFR.np

    CASE = 3;
    mesh_dim = size(varargin{1});

    samples = cellfun(@(s) {s(:)'}, varargin);
    samples = vertcat(samples{:});

end

if CASE == 0
    error('Arguments are not good: numel(varargin) = %d, np = %d.', numel(varargin), pLFR.np)
end

% samples = [
%     1 2 3 4 1 2 3
%     [ 1 2 3 4 1 2 3 ] * 2
%     [ 1 2 3 4 1 2 3 ] * 3
%     ];

%  [ { p11   p21   ...  pN1  }
%    { ...   ...   ...  ...  }
%    { p1np  p2np  ...  pNnp } ] --> celling along 2nd dim
%

samples_cell = num2cell(samples,2);

% Egy nullaval kiegeszitve: 2019.11.06. (november  6, szerda), 01:39
delta = pLFR.delta_fh(samples_cell{:},zeros(size(samples_cell{1})));

delta_cell = num2cell(delta,1);

val = cellfun( @(d) { pLFR.A + pLFR.B/( pLFR.I - diag(d)*pLFR.D)*diag(d)*pLFR.C }, delta_cell);


switch CASE
    case 1
        varargout{1} = horzcat(val{:});
        
    case 2
        varargout{1} = vertcat(val{:});
        
    case 3
        varargout{1} = reshape(vertcat(val{:}),mesh_dim);
        
end


end