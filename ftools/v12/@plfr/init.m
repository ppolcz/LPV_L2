function pLFR = init(pLFR,varargin) %(arg1,arg2,arg3,arg4,arg5,arg6)
%%

NrIn = nargin-1;

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
% Created:    2019.11.20. (november 20, szerda), 15:30
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

    case 1 % (syslfr) => generate_Delta
        pLFR.lfrtbx_obj = varargin{1};
        [pLFR.D,pLFR.C,pLFR.B,pLFR.A] = lfrdata(pLFR.lfrtbx_obj);


    case 2 % (M,blk) => generate_Delta
        M = varargin{1};
        blk = varargin{2};
        m1 = sum(blk.desc(1,:));

        [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = pcz_split_matrix(M, [Inf m1], [Inf m1]);
        pLFR.lfrtbx_obj = lfr(pLFR.D,pLFR.C,pLFR.B,pLFR.A,blk);


    case 3 % (M,Delta,bounds) => generate_LFR
        M = varargin{1};
        pLFR.Delta = varargin{2};
        bounds = varargin{3}';
        m1 = size(pLFR.Delta,1);

        [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = pcz_split_matrix(M, [Inf m1], [Inf m1]);


    case 5 % (A,B,C,D,blk) => generate_Delta
        pLFR.lfrtbx_obj = lfr(varargin{[4,3,2,1,5]});
        [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = deal(varargin{1:4});


    case 6 % (A,B,C,D,5:Delta,6:bounds) => generate_LFR
        [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = deal(varargin{1:4});
        pLFR.Delta = varargin{5};
        bounds = varargin{6}';
        
    % This case is more special compared to the others.
    case 7 % (A,B,C,D,5:Delta,6:[],7:blk) => generate_LFR
        [pLFR.A,pLFR.B,pLFR.C,pLFR.D] = deal(varargin{1:4});
        pLFR.Delta = varargin{5};
        pLFR = pLFR.generate_LFR(varargin{7});

end

% Generate (Delta and bounds) or lfrtbx_obj
switch NrIn

    case {1,2,5} % lfrtbx_obj --> (Delta, bounds)
        pLFR = pLFR.generate_Delta;

    case {3,6} % (Delta, bounds) --> (lfrtbx_obj, bounds -- update)
        pLFR = pLFR.generate_LFR(bounds);

end


% Modified [BEGIN], 2020.03.29. (március 29, vasárnap), 15:50

pLFR = pLFR...
    .generate_dimensions...
    .generate_other_fields;

if isempty(req_subsvars)
    pLFR = pLFR.sort_vars;
else
    pLFR = pLFR.set_vars(req_subsvars(:));
end

pLFR = pLFR.test_properties;

% Modified [END]

end

