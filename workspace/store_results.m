function store_results(fname, modelname, gamma_lower, gamma_upper, solver_time, overall_time, info, method)
%  
%  File: store_results.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2020. March 25. (2019b)

RUN_ID = str2double(getenv('RUN_ID'));
if ~RUN_ID
    RUN_ID = pcz_runID;
end

s.DateTime = datestr(now, 'yyyy.mm.dd. dddd HH:MM:SS');
s.Model = modelname;
s.RunID = RUN_ID;
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

[~,fname,~] = fileparts(fname);

Results_spreadsheet = [ 'results' filesep fname '.xlsx' ];

if ~exist('results','dir')
    mkdir result
end

if exist(Results_spreadsheet, 'file')
    Results = readtable(Results_spreadsheet,'Sheet',1);
end

if ~exist('Results', 'var')
    Cols = { 'DateTime', 'Model', 'RunID', 'Lower', 'Upper', 'Solver_Time', 'Overall_Time', 'Solver_Info', 'Method' };
    Results = cell2table(cell(0,numel(Cols)), 'VariableNames', Cols);
end

Results = [ Results ; struct2table(s) ];

writetable(Results,Results_spreadsheet,'Sheet',1);

pcz_dispFunction('Results stored in `%s''', Results_spreadsheet);
pcz_dispFunction2(evalc('disp(s)'))


latex_fname = getenv('LATEX_FNAME');
if isempty(latex_fname)
    latex_fname = 'latex_source_code.tex';
end


latexfile = fopen(latex_fname,'a');

fprintf(latexfile,'\n%s %s \n\t& %.4g & %d & %d & %d \\\\\n',...
    method,info,gamma_upper,round(solver_time),round(overall_time),round(overall_time) - round(solver_time));

fclose(latexfile);

end