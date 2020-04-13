function helper_generate_isosurface_method5(Q,PI_x,gamma,p_expr,x_lim,p_lim,dp_lim,resolution)
%
%  File: helper_generate_isosurface_method5.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2020. March 26. (2019b)
%  Major review on 2020. April 09. (2019b)

%%

TMP_oyXrvtrEzBjWNTyEyMog = pcz_dispFunctionName('Generate iso-surface of the Lyapunov/storage function');

np = size(p_expr,1);
nx = size(PI_x,2);
P_generate_symvars(nx,np,0,0);

TMP_DGzdAcREJHGXjyRzNjJV = pcz_dispFunctionName('Conversion to `plfr'' object');

[x_lfr,x_lfr_cell] = pcz_generateLFRStateVector('x',x_lim);
[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim,dp_lim);


% pi(x,p) = Pi(p) * x as a `plfr' object
pi_x = PI_x * x_lfr;

% Q(p) as a `plfr' object
Qlfr = plfr(Q,p_lfr_cell);

% V(x,p) = pi'(x,p) Q(p) pi(x,p) as a `plfr' object
V_lfr = pi_x' * Qlfr * pi_x;

pcz_dispFunctionEnd(TMP_DGzdAcREJHGXjyRzNjJV);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

TMP_iwclITrbVyVrJaArrXNr = pcz_dispFunctionName('Evaluate the storage function over the grid');

% Generate grid in X
lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([x_lim resolution(:)]),2));
x_grid = cell(1,nx);
[x_grid{:}] = ndgrid(lspace{:});

% Vectorize grid in X
grid_size = size(x_grid{1});
xx = cellfun(@(a) {a(:)}, x_grid);

% Compute the values of p over the grid in X
p_fh = matlabFunction(p_expr(:).','vars',x);
pp = p_fh(xx{:});

% Value of (x,p) over the vectorized grid in X
xxpp = horzcat(xx{:},pp);

% Compute the values of the storage function
VV = V_lfr.val(xxpp);
V_num = reshape(VV,grid_size);

pcz_dispFunctionEnd(TMP_iwclITrbVyVrJaArrXNr);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

TMP_GvDXGhRLfipwBoRPoGfI = pcz_dispFunctionName('Find maximal level set inside X');

% Find maximal level set inside X
alpha = sdpvar;

% Hard coded (Now, I am a too lazy, to find a more adequate solution)
CONS = [
    V_num([1 end],:,:) >= alpha+eps
    V_num(:,[1 end],:) >= alpha+eps
    V_num(:,:,[1 end]) >= alpha+eps
    ];

sol_alpha = optimize(CONS,-alpha,sdpsettings('verbose',0));
alpha = double(alpha);

Mu = alpha/gamma;

pcz_dispFunction2('alpha = %g', alpha)
pcz_dispFunction2('gamma = %g', gamma)
pcz_info('|u|_L2 < M = alpha/gamma = %g', Mu)

pcz_2basews(alpha,Mu)

pcz_dispFunctionEnd(TMP_GvDXGhRLfipwBoRPoGfI);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if alpha > eps
%% PLOT if a good level set obtained

    figure('Color', [1 1 1])
    ph = patch(isosurface(x_grid{:},V_num,alpha));
    ph_Vertices = ph.Vertices;
    ph_Faces = ph.Faces;
    set(ph,'FaceColor','red','EdgeColor','none', 'FaceAlpha', 0.8);


    axis equal
    light('Position',[-5 0 0]),
    view([-14 15])
    xlabel $x_1$ interpreter latex
    ylabel $x_2$ interpreter latex
    zlabel $x_3$ interpreter latex

    grid on
    axlims = x_lim';
    axis(axlims(:)),
    box on
    % set(gca,'LineWidth',1)

    set(gca, Logger.font_axis12{:})
    Logger.latexify_axis(gca,20)
    Labels = Logger.latexified_labels(gca,28,'$x_1 = v$','$x_2 = \vartheta$','$x_3 = \omega$');
    Labels{1}.VerticalAlignment = 'middle';
    Labels{2}.VerticalAlignment = 'bottom';

    good = V_num < alpha;
    all = ones(size(V_num));

    VOLUME = prod(x_lim(:,2) - x_lim(:,1)) * ( sum(good(:))/sum(all(:)) );

	pcz_dispFunction2('Approximated volume of the invariant domain: %g', VOLUME)

    rotate3d on

else

    pcz_info(0,'Iso surface completely inside X not found')

end

pcz_dispFunctionEnd(TMP_oyXrvtrEzBjWNTyEyMog);

end
