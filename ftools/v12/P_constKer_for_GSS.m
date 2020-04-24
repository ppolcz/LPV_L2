function [constKer,M,samples,sampleserr,error_] = P_constKer_for_GSS(M_gss,varargin)
%% P_constKer_for_GSS
%
%  File: P_constKer_for_GSS.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. March 10.
%  Major review on 2020. March 11. (2019b, v11)
%  Major review on 2020. April 10. (2019b, v12)
%

%%

if isa(M_gss,'plfr')
    M_gss = M_gss.lfrtbx_obj;
end

if isa(M_gss,'lfr')
    try
        M_gss = lfr2gss(M_gss);
    catch e
        if size(M_gss,'order') == 0
            M_gss = gss(M_gss.d);
        else
            error(e)
        end
    end          
end

args.lims = []; % make it back-compatible
args.tol = 1e-7;
args.proj = [];
args.samples = [];
args.sampleserr = [];
args.N = 0;
args = parsepropval(args,varargin{:});

DO_PROJECTION = true;
if isempty(args.proj)
    DO_PROJECTION = false;
    args.proj = @(W) W;
end

if DO_PROJECTION
    error('P_constKer_for_GSS cannot handle projection yet, try P_constKer_for_LFR instead!')
end

%%

TMP_RNSQSkSbcBzsqumHhrSV = pcz_dispFunctionName;

[ny,nu] = size(M_gss);

pcz_dispFunction2('Constant kernel computation for G: %dx%d, with tolerance: %g.',...
    ny,nu,args.tol)
pcz_dispFunctionSeparator;

Optimal_N = ceil(nu / ny);
N = max([args.N,Optimal_N,2]);

M = zeros(0,nu);

if ~isempty(args.samples) && ~isempty(args.sampleserr)
    %% If the samples are given preliminarily

    warning 'NOT TESTED: 2020.04.10. (április 10, péntek), 20:50'
    
    .......................................................................
    ... MAIN COMPUTATIONS .................................................
    .......................................................................

    % Load samples for kernel computation if it is given.
    samples = args.samples;
    Nr_Samples = size(samples{2},1);

    M = eval(M_gss,samples{:});
    M = vertcat(M{:});
    constKer = null(M);

    pcz_dispFunction2('Kernel computed for M: %dx%d, where (%d = %d*%d)',...
        Nr_Samples*ny,nu,Nr_Samples*ny,Nr_Samples,ny);


    .......................................................................
    ... ERROR COMPUTATIONS ................................................
    .......................................................................

    % Load samples for error computation if it is given.
    sampleserr = args.sampleserr;
    N = size(sampleserr{2},1);

    Merr = eval(M_gss, sampleserr{:});

    M = [ M ; vertcat(Merr{:})];
    Nr_Samples = Nr_Samples + N;

    error_ = norm(M * constKer);

    pcz_dispFunction2('Error computed for M: %dx%d, where (%d = %d*%d). Error: %d',...
        Nr_Samples*ny,nu,Nr_Samples*ny,Nr_Samples,ny, error_);


    .......................................................................
    ... CHECK SOLUTION ....................................................
    .......................................................................


    % norm_LFR0 = distlfr(LFR0, minlfr(LFR0*0));
    if error_ < args.tol
        return;
    end

    pcz_dispFunction2('Constant kernel with N = %d not found, trying again with N = %d.', ...
        Nr_Samples-N, Nr_Samples+N)
    pcz_dispFunction2('Though the samples were given by the user.')

end

sampleserr = [];
error_ = [];
Nr_Samples = N;
[M,samples] = dbsample(M_gss,N);
M = vertcat(M{:});
constKer = null(M);

Nr_iterations = ceil(10 * Optimal_N / N);
i = 1;
while true
    %% If the samples are not given, nor the limits: generate random points

    pcz_dispFunction2('Kernel computed for M: (%d*%d=)%dx%d',...
        Nr_Samples,ny,size(M));

    if isempty(constKer)
        pcz_info(1,'Kernel = {0}.')
        break
    end
    .......................................................................
    ... ERROR COMPUTATION .................................................
    .......................................................................

    ZERO_gss = mingss(M_gss*constKer);
    
    [~,~,m1,~,~] = size(ZERO_gss);
    if m1 == 0
        ZERO = ZERO_gss.M.D;
    else
        [ZERO_sampled,~] = dbsample(ZERO_gss,N);
        ZERO = vertcat(ZERO_sampled{:});        
    end

    % Taking the maximum value
    error_ = max(abs(ZERO(:)));

    pcz_dispFunction2('Error computed for Merr: (%d*%d=)%dx%d. Maximum error: %d',...
        N,ny,N*ny,nu,error_);

    
    .......................................................................
    ... DISPLAY STATUS ....................................................
    .......................................................................

    pcz_dispFunctionSeparator;

    if error_ < args.tol
        pcz_info(true,'Kernel found with the given tolerance.')
        break;
    elseif i >= Nr_iterations
        pcz_info(false,'Maximum nr. iterations (%d) reached. Kernel NOT found with the given tolerance.', ...
            Nr_iterations)
        break;
    else
        i = i + 1;
    end
    
    .......................................................................
    ... MAIN COMPUTATIONS .................................................
    .......................................................................


    pcz_dispFunction2('Constant kernel with N = %d not found, trying again with N = %d.', ...
        Nr_Samples-N, Nr_Samples+N)

    [M_sampled,new_samples] = dbsample(M_gss,N);
    samples = [samples new_samples];

    M = [ M ; vertcat(M_sampled{:}) ];
    Nr_Samples = Nr_Samples + numel(new_samples);

    constKer = null(M);

end

pcz_dispFunctionEnd(TMP_RNSQSkSbcBzsqumHhrSV);

end

function test1_simple_example
%%
    G_reset

    [x,x_cell] = pcz_generateGSStateVector('x',2)

    PI_fh = @(x1,x2) [
        1 0
        0 1
        x1 0
        0 x2
        x1^2*x2 0
        0 x1*x2
        x1*x2^2 x1*x2
        ];

    PI = PI_fh(x_cell{:});

    [m,n] = size(PI);

    N_Xi = [
        eye(m)
        x_cell{1}*eye(m)
        x_cell{2}*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0,M] = P_constKer_for_GSS(Aminek_keresem','N',1);

    [U,Sigma,V] = svd(M);

    Size_N0 = size(N0)
    Size_M = size(M)

    format long g
    diag(Sigma)
    format short

    ZERO = N0' * Aminek_keresem;

    pcz_gsszero_report(ZERO)

end

function test2_trivial_kernel
%%

    G_reset

    TMP_LFTwIriJzvmzScZYhZgI = pcz_dispFunctionName('test2_trivial_kernel');

    [~,x_cell] = pcz_generateGSStateVector('x',3);

    PI_fh = @(x1,x2,x3) [
        1 0 0
        0 1 0
        x1 0 1
        0 x2 x1*x3
        x1^2*x2 0 x3
        0 x1*x2 x2*x3^2
        x1*x2^2 x1*x2 1/(x1^2 + x3^2 + 1)
        ];
    PI = PI_fh(x_cell{:});

    [m,n] = size(PI);

    % N_XI = kronLFR([1;x1;x2],eye(m));
    N_Xi = [
        eye(m)
        x_cell{1}*eye(m)
        x_cell{2}*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0,M] = P_constKer_for_GSS(Aminek_keresem');

    [U,Sigma,V] = svd(M);

    Size_N0 = size(N0)
    Size_M = size(M)

    format long g
    diag(Sigma)
    format short

    ZERO = N0' * Aminek_keresem;

    pcz_gsszero_report(ZERO)

    pcz_dispFunctionEnd(TMP_LFTwIriJzvmzScZYhZgI);

end

function test3
%%

G_reset



P_generate_symvars_v10(8,0,0,0);

x_lim = repmat([1,3],[8 1]);

[~,x_cell] = pcz_generateGSStateVector('x',x_lim);

PI_fh = @(x1,x2,x3,x4,x5,x6,x7,x8) [
    1 0 0
    0 x8 0
    x5*x3^2/(1 + 72^2*x1^2 + x8*x1) 0 x1*x2
    1/(x1^2 + x3^2 + 1) x1*x2 1
    x1 0 1
    0 x2 x1*x3
    x1^2*x2*(1+x2^2*x3)/(6+x7) 0 x3
    0 x1*x5 x2*x3^2/(1 + x4^2*x1^2 + x2*x1)
    x1*x2^2 x1*x2 1/(x1^2 + x3^2 + 1)
    ];
PI = PI_fh(x_cell{:});

[m,n] = size(PI);

S = randn(m+4,m);

PI = S*PI;
[m,n] = size(PI);

% N_XI = kronLFR([1;x1;x2],eye(m));
N_Xi = [
    eye(m)
    x_cell{1}*eye(m)
    x_cell{2}*eye(m)
    x_cell{3}*eye(m)
    x_cell{4}*eye(m)
    x_cell{5}*eye(m)
    x_cell{6}*eye(m)
    x_cell{7}*eye(m)
    x_cell{8}*eye(m)
    x_cell{1}^2*eye(m)
    x_cell{1}*x_cell{2}*eye(m)
    x_cell{2}^2*eye(m)
    x_cell{1}*x_cell{3}*eye(m)
    x_cell{2}*x_cell{3}*eye(m)
    x_cell{3}^2*eye(m)
    x_cell{5}^2*eye(m)
    x_cell{5}*x_cell{6}*eye(m)
    x_cell{6}^2*eye(m)
    x_cell{5}*x_cell{7}*eye(m)
    x_cell{6}*x_cell{7}*eye(m)
    x_cell{7}^2*eye(m)
    ];

Aminek_keresem = N_Xi * PI;

[N0_v1,M] = P_constKer_for_GSS(Aminek_keresem');

[N0_v2,M] = P_constKer_for_GSS(Aminek_keresem','N',20);

[U,Sigma,V] = svd(M);

Size_M = size(M)

format long g
diag(Sigma)
format short

ZERO = N0_v1' * Aminek_keresem;
pcz_gsszero_report(ZERO,x,'Elso')

ZERO = N0_v2' * Aminek_keresem;
pcz_gsszero_report(ZERO,x,'Masodik')

Size_N0_v1 = size(N0_v1)
Size_N0_v2 = size(N0_v2)

end
