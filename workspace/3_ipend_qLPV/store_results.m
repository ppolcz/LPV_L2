function store_results(fname, x_lim, gamma_lower, gamma_upper, solution_time, info, method, Use_MinLFR)
%  
%  File: store_results.m
%  Directory: 1_PhD_projects/22_Hinf_norm/qLPV_ipend
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 14. (2019b)
%

s.Stamp = datestr(now, 'yyyy.mm.dd. dddd HH:MM:SS');
s.x1_max = x_lim(1,2);
s.x2_max = x_lim(2,2);
s.x3_max = x_lim(3,2);
s.Lower = gamma_lower;
s.Upper = gamma_upper;
s.Solver_Time = solution_time;
s.Solver_info = info;
s.Method = method;
s.Use_MinLFR = Use_MinLFR;

if isempty(s.Upper)
    s.Upper = 0;
end

Results_csv = [ 'results/' fname];

if exist(Results_csv, 'file')
    Results = readtable(Results_csv,'Delimiter','|');
end

if ~exist('Results', 'var')
    Cols = { 'Stamp', 'x1_max', 'x2_max', 'x3_max', 'Lower', 'Upper', 'Solver_Time', 'Solver_info', 'Method', 'Use_MinLFR' };
    Results = cell2table(cell(0,numel(Cols)), 'VariableNames', Cols);
end

Results = [ Results ; struct2table(s) ];

% save(Results_mat,'Results');
writetable(Results,Results_csv,'Delimiter','|');

end