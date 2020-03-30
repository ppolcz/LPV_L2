function store_results(fname, modelname, gamma_lower, gamma_upper, solver_time, overall_time, info, method)
%  
%  File: store_results.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2020. March 25. (2019b)

RUN_ID = getenv('RUN_ID');
if ~RUN_ID
    RUN_ID = num2str(pcz_runID);
end

s.Model = modelname;
s.Stamp = [ datestr(now, 'yyyy.mm.dd. dddd HH:MM:SS') ' id:' RUN_ID ];
s.Lower = gamma_lower;
s.Upper = gamma_upper;
s.Solver_Time = solver_time;
s.Overall_Time = overall_time;
s.Solver_Info = info;
s.Method = method;

if isempty(s.Lower)
    s.Lower = 0;
end

if isempty(s.Upper)
    s.Upper = 0;
end

Results_csv = [ 'results/' fname];

if exist(Results_csv, 'file')
    Results = readtable(Results_csv,'Delimiter','|');
end

if ~exist('Results', 'var')
    Cols = { 'Model', 'Stamp', 'Lower', 'Upper', 'Solver_Time', 'Overall_Time', 'Solver_Info', 'Method' };
    Results = cell2table(cell(0,numel(Cols)), 'VariableNames', Cols);
end

Results = [ Results ; struct2table(s) ];

writetable(Results,Results_csv,'Delimiter','|');

pcz_dispFunction('Results stored in `%s''', Results_csv);
pcz_dispFunction2(evalc('disp(s)'))

end