function helper_generate_isosurface_method5(Q,PI_x,gamma,p_expr,x_lim,resolution)
%
%  File: helper_generate_isosurface_method5.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2020. March 26. (2019b)

%%
TMP_oyXrvtrEzBjWNTyEyMog = pcz_dispFunctionName('Generate iso-surface of the Lyapunov/storage function');

np = size(p_expr,1);
nx = size(PI_x,2);
P_generate_symvars(nx,np,0,0);

Qx = Q.set_channels([1 p_expr.'],x);

PI_x_subs = subs(PI_x(p_expr));

V = matlabFunction(x.' * PI_x_subs.' * sym(Qx) * PI_x_subs * x);

resolution = resolution(:);

lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([x_lim resolution]),2));
xx = cell(1,nx);
[xx{:}] = ndgrid(lspace{:});

VV = V(xx{:});

alpha = sdpvar;

% Hard coded (Now, I am a too lazy, to find a more adequate solution)
CONS = [
    VV([1 end],:,:) >= alpha+eps
    VV(:,[1 end],:) >= alpha+eps
    VV(:,:,[1 end]) >= alpha+eps
    ];

sol_alpha = optimize(CONS,-alpha,sdpsettings('verbose',0));
alpha = double(alpha);

Mu = alpha/gamma;

pcz_dispFunction2('alpha = %g', alpha)
pcz_dispFunction2('gamma = %g', gamma)
pcz_info('|u|_L2 < M = alpha/gamma = %g', Mu)

pcz_2basews(alpha,Mu)

if alpha > eps

    figure('Color', [1 1 1])
    ph = patch(isosurface(xx{:},VV,alpha));
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

    good = VV < alpha;
    all = ones(size(VV));

    VOLUME = prod(x_lim(:,2) - x_lim(:,1)) * ( sum(good(:))/sum(all(:)) );

	pcz_dispFunction2('Approximated volume of the invariant domain: %g', VOLUME)

    rotate3d on

else

    pcz_info(0,'Iso surface completely inside X not found')

end

pcz_dispFunctionEnd(TMP_oyXrvtrEzBjWNTyEyMog);

end
