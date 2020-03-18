function [v, ak, fci, proj_x, proj_xp, ek, Vk] = P_ndnorms_of_X(x_lim, p_lim, varargin)
%% 
%  
%  file:   P_ndnorms_of_X.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.04.17. Sunday, 16:40:39
%  Modified on 2019. November 19. (2019a) -- projection extended xp space
%  Modified on 2020. January 05. (2019b) -- additional segments (multiple)
%

% 2020.01.05. (január  5, vasárnap), 00:32 [INNEN]
np = 0;
if nargin >= 2 && ischar(p_lim)
    varargin = [ p_lim varargin ];
elseif nargin >= 2
    np = size(p_lim,1);
end

args.Vertices = [];
args.Facets = [];
args = parsepropval(args,varargin{:});
% 2020.01.05. (január  5, vasárnap), 00:32 [IDAIG]


% 2020.01.05. (január  5, vasárnap), 01:13
if ~isempty(x_lim)
    nx = size(x_lim,1);

    %% [1] Generate v
    x_lim_cell = num2cell(x_lim,2);
    v = allcomb(x_lim_cell{:});

    x_lim_bool = num2cell(repmat([0 1], nx, 1), 2);
    v_bool = allcomb(x_lim_bool{:});

    %% [2] Generate fci

    Nr_vertices_per_facet = 2^(nx-1);

    fci = zeros(nx*2,Nr_vertices_per_facet);

    for i = 1:nx
        for k = [0 1]
            fci(2*i+k-1,:) = find(v_bool(:,i) == k)';
        end
    end
else
    v = [];
    fci = [];
end

% figure,
% trisurf(fci(:,[1 2 4 3]),v(:,1), v(:,2), v(:,3), 'Facecolor','blue', 'facealpha', 0.5); axis equal;

%% [1.1] [2.1] Kiegeszites

% 2020.01.15. (január 15, szerda), 19:27 [INNEN]
if ~isempty(args.Vertices) && size(args.Vertices,2) == 2 && isempty(args.Facets)
    
    I = convhull(args.Vertices);
    args.Facets = [
        I(1:end-1) I(2:end)
        ];
    
end
% 2020.01.15. (január 15, szerda), 19:27 [IDAIG]

% 2020.01.05. (január  5, vasárnap), 00:32 [INNEN]
if ~isempty(args.Vertices) && ~isempty(args.Facets)

    fci = [
        fci
        args.Facets + size(v,1)
        ];
    
    v = [
        v
        args.Vertices
        ];

end
% 2020.01.05. (január  5, vasárnap), 00:32 [IDAIG]

%% [3] Generate ek, ak and Vk

nx = size(v,2);
Nr = size(fci,1);
ek = zeros(Nr,nx);
ak = zeros(Nr,nx);
Vk = cell(Nr,1);
proj_x = cell(Nr,1);
proj_xp = cell(Nr,1);

Nr_vertices_per_facet = 2^(nx-1);

for i = 1:Nr

    r0 = v(fci(i,1), :)';

    % 2D: W = [ B - A ] - actually calculated
    % 3D: W = [ B - A ; C - A ; D - A ] - actually calculated
    % 3D: W = [ B - A ; D - A ] - it would be enough
    % 4D: W = [ A2 - A1 ; B1 - A1 ; B2 - A1 ; C1 - A1 ; C2 - A1 ; D1 - A1 ; D2 - A1 ] - actually calculated
    % 4D: W = [ A2 - A1 ; B1 - A1 ; C1 - A1 ; D1 - A1 ] - it would be enough
    V = (  v(fci(i,2:Nr_vertices_per_facet),:) - repmat(r0', [Nr_vertices_per_facet-1, 1])  )';
    
    % normal vector (norm(n) == 1)
    nvec = null(V');

    % [3.1] ek
    ek(i,:) = nvec;

    % distance of the origin to the dim dimensional facet
    dst = dot(nvec, r0);
    
    % [3.2] ak (may be undefined - division by 0)
    ak(i,:) = nvec / dst;
    
    % [3.3] Vk
    V = orth(V);
    Vk{i} = V;
    
    % [3.4.1] projection (only x needed to be substituted)
    proj_x{i} = @(x) V*V'*x + dst*nvec*ones(1,size(x,2));
    
    V_xp = blkdiag(V,eye(np));
    nvec_xp = [nvec ; zeros(np,1)];

    % [3.4.2] projection in extended space ([x;p] needed to be substituted)
    proj_xp{i} = @(xp) V_xp*V_xp'*xp + dst*nvec_xp*ones(1,size(xp,2));
end


end


function test1
    %%
    
    x_lim = [
        -2 1 
        -4 3
        ];
    
    [v,ak,fci,proj] = P_ndnorms_of_X(x_lim)
    
    
    for i = 1:size(fci,1)
    
        W = randn(2,50);
        Wi = proj{i}(W);
    
        pcz_plotpoint(Wi), hold on;
        
    end
    
    axis equal
    axis tight
    
end

function test2
    %%
    
    x_lim = [
        -1 1 
        -1 1
        -1 1
        ];
    
    [v,ak,fci,proj] = P_ndnorms_of_X(x_lim)
    
    
    for i = 1:size(fci,1)
    
        W = randn(3,50);
        Wi = proj{i}(W);
    
        pcz_plotpoint(Wi), hold on;
        
    end
    
    axis equal
    axis tight
    
end