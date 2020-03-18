function [constKer,M,samples,sampleserr,error_LFR0] = P_constKer_for_LFR(pLFR,varargin)
%% P_constKer_for_LFR
%  
%  File: P_constKer_for_LFR.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. March 10.
%  Major review on 2020. March 11. (2019b)
%

%%

if isa(pLFR,'lfr')
    pLFR = plfr(pLFR);
end
np = pLFR.np;

args.proj = @(W) W;
args.N = 0;
args.samples = [];
args.sampleserr = [];
args.lims = zeros(np,0);
args = parsepropval(args,varargin{:});

%%

TMP_RNSQSkSbcBzsqumHhrSV = pcz_dispFunctionName;

Optimal_N = ceil(pLFR.nu / pLFR.ny);

if args.N == 0
    N = Optimal_N;
else
    N = args.N;
end

Nr_Samples = 0;
M = zeros(0,pLFR.nu);

samples = zeros(np,0);
sampleserr = zeros(np,0);

if ~isempty(args.samples) && ~isempty(args.sampleserr)
    %% If the samples are given preliminarily

    .......................................................................
    ... MAIN COMPUTATIONS .................................................
    .......................................................................
    
    % Load samples for kernel computation if it is given.
    samples = args.samples;
    Nr_Samples = size(samples,2);

    M = pLFR.val(samples');
    constKer = null(M);

    pcz_dispFunction('Kernel computed for M: %dx%d, where (%d = %d*%d)',...
        Nr_Samples*pLFR.ny,pLFR.nu,Nr_Samples*pLFR.ny,Nr_Samples,pLFR.ny);
    
    
    .......................................................................
    ... ERROR COMPUTATIONS ................................................
    .......................................................................
    
    % Load samples for error computation if it is given.
    sampleserr = args.sampleserr;
    N = size(sampleserr,2);

    % 2019.09.05. (szeptember  5, csütörtök), 16:17
    M = [ M ; pLFR.val(sampleserr')];
    Nr_Samples = Nr_Samples + N;

    error_LFR0 = norm(M * constKer);
    
    pcz_dispFunction('Error computed for M: %dx%d, where (%d = %d*%d). Error: %d',...
        Nr_Samples*pLFR.ny,pLFR.nu,Nr_Samples*pLFR.ny,Nr_Samples,pLFR.ny, error_LFR0);

    
    .......................................................................
    ... CHECK SOLUTION ....................................................
    .......................................................................
    
    
    % norm_LFR0 = distlfr(LFR0, minlfr(LFR0*0));
    if error_LFR0 < 1e-7
        return;
    end

    pcz_dispFunction('Constant kernel with N = %d not found, trying again with N = %d.', ...
        Nr_Samples-N, Nr_Samples+N)
    pcz_dispFunction('Though the samples were given by the user.')
    
end


Nr_iterations = ceil(10 * Optimal_N / N);
for i = 1:Nr_iterations
    %% If the samples are not given, nor the limits: generate random points

    .......................................................................
    ... MAIN COMPUTATIONS .................................................
    .......................................................................
    
    % Generate new samples (if not given)
    new_samples = generate(N,args.lims,args.proj);
    samples = [ samples sampleserr new_samples ];
    

    M = [ M ; pLFR.val(new_samples')];
    Nr_Samples = Nr_Samples + size(new_samples,2);

    constKer = null(M);

    pcz_dispFunction('Kernel computed for M: %dx%d, where (%d = %d*%d)',...
        size(M),Nr_Samples*pLFR.ny,Nr_Samples,pLFR.ny);
    
    
    .......................................................................
    ... ERROR COMPUTATIONS ................................................
    .......................................................................
            
    sampleserr = generate(N,args.lims,args.proj);
            
    M = [ M ; pLFR.val(sampleserr')];
    Nr_Samples = Nr_Samples + size(sampleserr,2);

    error_LFR0 = norm(M * constKer);
    
    pcz_dispFunction('Error computed for M: %dx%d, where (%d = %d*%d). Error: %d',...
        Nr_Samples*pLFR.ny,pLFR.nu,Nr_Samples*pLFR.ny,Nr_Samples,pLFR.ny, error_LFR0);

    
    .......................................................................
    ... CHECK SOLUTION ....................................................
    .......................................................................
    
    
    LFR0 = pLFR * constKer;
    LFR0 = minlfr(LFR0.lfrtbx_obj);
    if isempty(LFR0.d)
        break;
    end

    % norm_LFR0 = distlfr(LFR0, minlfr(LFR0*0));
    if error_LFR0 < 1e-7
        break;
    end

    pcz_dispFunction('Constant kernel with N = %d not found, trying again with N = %d.', ...
        Nr_Samples-N, Nr_Samples+N)
        
    pcz_dispFunctionSeparator;
    
end

pcz_dispFunctionEnd(TMP_RNSQSkSbcBzsqumHhrSV);

end

function samples = generate(N,lims,proj)

    if isempty(lims)
        samples = proj(randn(size(lims,1),N));
    else    
        samples_cell = cellfun(@(a,b) {rand(1,N)*(b-a)+a}, num2cell(lims(:,1)), num2cell(lims(:,2)));
        samples = proj(vertcat(samples_cell{:}));
        
        % After the projection some points may lie outside `lims'.
        samples = samples(:,prod(lims(:,1)*ones(1,N) <= samples & samples <= lims(:,2)*ones(1,N),1) == 1);
        
        %{
        plot3(samples(1,:),samples(2,:),samples(3,:),'.','MarkerSize',10)
        %}
        
    end

end


function test1_simple_example
%%

    global SCOPE_DEPTH VERBOSE LATEX_EQNR 
    SCOPE_DEPTH = 0;
    VERBOSE = 1;
    LATEX_EQNR = 0;

    x = sym('x',[2 1]);

    lfrs x1 x2 
    PI = [
        1 0 
        0 1 
        x1 0
        0 x2
        x1^2*x2 0
        0 x1*x2
        x1*x2^2 x1*x2
        ];

    [m,n] = size(PI);

    % N_XI = kronLFR([1;x1;x2],eye(m));
    N_Xi = [
        eye(m)
        x1*eye(m)
        x2*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0,~,M] = P_constKer_for_LFR(Aminek_keresem','N',1);

    [U,Sigma,V] = svd(M);

    Size_N0 = size(N0)
    Size_M = size(M)

    format long g
    diag(Sigma)
    format short
    
    ZERO = plfr(N0' * Aminek_keresem);
    
    pcz_fhzero_report(ZERO,x)

end

function test2_trivial_kernel
%%

    x = sym('x',[3 1]);

    lfrs x1 x2 x3
    PI = [
        1 0 0
        0 1 0
        x1 0 1
        0 x2 x1*x3
        x1^2*x2 0 x3
        0 x1*x2 x2*x3^2
        x1*x2^2 x1*x2 1/(x1^2 + x3^2 + 1)
        ];

    [m,n] = size(PI);

    % N_XI = kronLFR([1;x1;x2],eye(m));
    N_Xi = [
        eye(m)
        x1*eye(m)
        x2*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0,~,M] = P_constKer_for_LFR(Aminek_keresem');

    [U,Sigma,V] = svd(M);

    Size_N0 = size(N0)
    Size_M = size(M)

    format long g
    diag(Sigma)
    format short

    ZERO = plfr(N0' * Aminek_keresem);
    
    pcz_fhzero_report(ZERO,x)
    
end

function test3
%%

    global SCOPE_DEPTH VERBOSE LATEX_EQNR 
    SCOPE_DEPTH = 0;
    VERBOSE = 1;
    LATEX_EQNR = 0;

    x = sym('x',[8 1]);

    lfrs x1 x2 x3 x4 x5 x6 x7 x8
    PI = [
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
    [m,n] = size(PI);
    
    S = randn(m+4,m);
    
    PI = S*PI;
    [m,n] = size(PI);

    % N_XI = kronLFR([1;x1;x2],eye(m));
    N_Xi = [
        eye(m)
        x1*eye(m)
        x2*eye(m)
        x3*eye(m)
        x4*eye(m)
        x5*eye(m)
        x6*eye(m)
        x7*eye(m)
        x8*eye(m)
        x1^2*eye(m)
        x1*x2*eye(m)
        x2^2*eye(m)
        x1*x3*eye(m)
        x2*x3*eye(m)
        x3^2*eye(m)
        x5^2*eye(m)
        x5*x6*eye(m)
        x6^2*eye(m)
        x5*x7*eye(m)
        x6*x7*eye(m)
        x7^2*eye(m)
        ];

    Aminek_keresem = N_Xi * PI;

    [N0_v1,~,M] = P_constKer_for_LFR(Aminek_keresem');
    
    [N0_v2,~,M] = P_constKer_for_LFR(Aminek_keresem','N',20);

    [U,Sigma,V] = svd(M);

    Size_M = size(M)

    format long g
    diag(Sigma)
    format short

    ZERO = plfr(N0_v1' * Aminek_keresem);
    pcz_fhzero_report(ZERO,x,'Elso')

    ZERO = plfr(N0_v2' * Aminek_keresem);
    pcz_fhzero_report(ZERO,x,'Masodik')

    Size_N0_v1 = size(N0_v1)
    Size_N0_v2 = size(N0_v2)

end

% % 2020.03.11. (március 11, szerda), 23:13
% if ~isempty(args.lims)
%     %% Limits are given: use a grid! NOT GOOD!
%     
%     N = max([3^np,N])
%     
%     Ndim = ceil(N^(1/np));
%         
%     lspace_cell = cellfun(@(args) {linspace(args{:}, Ndim)}, num2cell(num2cell(args.lims),2));
%     samples_grid = cell(1,np);
%     [samples_grid{:}] = ndgrid(lspace_cell{:});
%     samples_cell = cellfun(@(s) {s(:)'}, samples_grid);
%     samples = args.proj(vertcat(samples_cell{:}));
%     
%     diff_dim = diff(args.lims') / Ndim;
%     ind_intm = cellfun(@(a) {1:Ndim-1}, num2cell(1:np));
%     sampleserr_grid = cellfun(@(s,offset) {s(ind_intm{:})+offset}, samples_grid, num2cell(diff_dim/2));
%     samples_cell = cellfun(@(s) {s(:)'}, sampleserr_grid);
%     sampleserr = args.proj(vertcat(samples_cell{:}));
%     
%     %{ 
%     figure, hold on
%     plot3(samples(1,:),samples(2,:),samples(3,:),'.','MarkerSize',100)
%     plot3(sampleserr(1,:),sampleserr(2,:),sampleserr(3,:),'.','MarkerSize',50)
%     %}  
%     
%     
%     .......................................................................
%     ... MAIN COMPUTATIONS .................................................
%     .......................................................................
%     
%     Nr_Samples = size(samples,2);
%     M = pLFR.val(samples');
%     constKer = null(M);
% 
%     pcz_dispFunction('Kernel computed over a grid for M: %dx%d, where (%d = %d*%d)',...
%         Nr_Samples*pLFR.ny,pLFR.nu,Nr_Samples*pLFR.ny,Nr_Samples,pLFR.ny);
%     
%     
%     .......................................................................
%     ... ERROR COMPUTATIONS ................................................
%     .......................................................................
%     
%     Nr_Samples = Nr_Samples + size(sampleserr,2);
%     M = [ M ; pLFR.val(sampleserr')];
% 
%     error_LFR0 = norm(M * constKer);
%     
%     pcz_dispFunction('Error computed in the interior points, M: %dx%d, where (%d = %d*%d). Error: %d',...
%         Nr_Samples*pLFR.ny,pLFR.nu,Nr_Samples*pLFR.ny,Nr_Samples,pLFR.ny,error_LFR0);
% 
%     
%     .......................................................................
%     ... CHECK SOLUTION ....................................................
%     .......................................................................
%     
%     
%     % norm_LFR0 = distlfr(LFR0, minlfr(LFR0*0));
%     if error_LFR0 < 1e-7
%         return;
%     end
% 
%     pcz_dispFunction('Constant kernel with N = %d not found, trying again with N = %d.', ...
%         Nr_Samples-N, Nr_Samples+N)
%     pcz_dispFunction('Though the samples were given by a grid.')
%     
%     pcz_dispFunctionSeparator;
% 
% end
