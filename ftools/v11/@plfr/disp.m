function disp(G)
%% 
%  File: disp.m
%  Directory: 7_ftools/ftools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 20. (2019a)
%

%%

%{

  plfr with properties:

    lfrtbx_obj: [8×1 lfr]
      delta_fh: @(x1,x2,p1,ZERO)[p1;x1;x1;x1;x2;x2;x2]
             M: [15×8 double]
             A: [8×1 double]
             B: [8×7 double]
             C: [7×1 double]
             D: [7×7 double]
         Delta: [7×7 sym]
             I: [7×7 double]
            np: 3
            nu: 1
            ny: 8
            m1: 7
      subsvars: [3×1 sym]
        bounds: [3×2 double]
           blk: [1×1 struct]
         names: {'p1'  'x1'  'x2'}
          vars: [3×1 sym]
       symvars: [1×3 sym]
          desc: [11×3 double]

%}

    [Sym,~,~,s] = symr(G,[]);
    if ~isempty(Sym)
        fprintf('Symbolical value (2-norm err of M = %d):\n\n', s.err)
        disp(Sym);
    end

    stringify = @(s) cell2mat(join(cellfun(@char, num2cell(s), 'UniformOutput', 0),', '));

    % Az LFR Toolbox object kimenetet egy kicsit megfaragom
    lfrvars = evalc('display(G.lfrtbx_obj)');
    lfrvars_ind = strfind(lfrvars,'LFR-');
    lfrvars = lfrvars(lfrvars_ind(1):end-1);
    lfrvars_ind = strfind(lfrvars,newline);
    lfrvars = lfrvars(lfrvars_ind+1:end);
    
    if ~isempty(G.names)
    
        % Ezt varom el: '1 (2x2), p1 (3x3), p2(2x2)'
        names = cell2mat(join(...
            cellfun(@(v,i,j) {sprintf('%s (%dx%d)',v,i,j)}, ...
                G.names,num2cell(G.desc(1,:)), num2cell(G.desc(2,:))),...
            ', '));

        fprintf('PLFR object (%d x %d), uncertainty block: (%d x %d)\n', ...
            G.ny, G.nu, sum(G.desc(1,:)), sum(G.desc(2,:)))
        fprintf('Dimensions: np = %d, ny = %d, nu = %d, m1 = %d.\n', ...
            G.np, G.ny, G.nu, G.m1)

        fprintf('-- \n%s\n', lfrvars)

        fhstr = char(G.delta_fh);
        fhstr_ind = strfind(fhstr,',ZERO)');
        fhargs = [ fhstr(1:fhstr_ind(1)-1) ')' ];

        fprintf('Names, dimensions:    %s \n', names)
        fprintf('Blocks (vars):        %s (%s)\n', stringify(G.vars), class(G.vars))
        fprintf('Symvars:              %s (%s, np = %d)\n', stringify(G.symvars), class(G.symvars), G.np)
        fprintf('Substitute into:      %s (%s)\n', stringify(G.subsvars), class(G.subsvars))
        fprintf('Function handle args: %s             TODO: too many redundant fields\n', fhargs)
    
    else
        
        fprintf('PLFR object (%d x %d), uncertainty block: (%d x %d)\n', ...
            G.ny, G.nu, sum(G.desc(1,:)), sum(G.desc(2,:)))
        
    end
end