classdef Logger
%%
%  Author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2018.02.11. Sunday, 11:17:40
%  Major review on 2020. March 17. (2019b) -- simplified
%

%%
properties (GetAccess = public, SetAccess = private)
    diary_fname;        
    latex_fname;
    results_fname;
    stamp;
    date;
end

methods (Access = public)

    function p = Logger(fname, varargin)

        RUN_ID = getenv('RUN_ID');

        p.stamp = datestr(now, 'mmddHHMMSS');
        p.date = datestr(now, 'yyyy.mm.dd. dddd, HH:MM:SS');
        
        [~,basename,~] = fileparts(fname);
        
        p.diary_fname = [ cd filesep 'results' filesep basename '-' p.stamp '_id' RUN_ID '-output.txt' ];
        p.latex_fname = [ cd filesep 'results' filesep basename '-' p.stamp '_id' RUN_ID '-latex.tex' ];
        p.results_fname = [ cd filesep 'results' filesep basename '-' p.stamp '_id' RUN_ID '-results.xlsx' ];
        
        setenv('DIARY_FNAME',p.diary_fname);
        setenv('LATEX_FNAME',p.latex_fname);
        setenv('RESULTS_FNAME',p.results_fname);
        
        [dirname,~,~] = fileparts(p.diary_fname);
        
        if ~exist(dirname,'dir')
            mkdir(dirname)
        end
        
        % Start diary (output logging)
        diary(p.diary_fname);
        pcz_info('Output logging (with `diary''): %s', p.diary_fname);

    end

    function stoplog(p)
        diary off
        if exist(p.diary_fname,'file')
            pcz_output2log(p.diary_fname);
            pcz_dispFunction2(' ')
            pcz_info('Logfile `%s'' formatted!', p.diary_fname);
        end
    end

end


properties (Constant = true)
    font_axis08 = {'FontSize',8,'FontName','TeX Gyre Schola Math'};
    font_axis10 = {'FontSize',10,'FontName','TeX Gyre Schola Math'};
    font_axis12 = {'FontSize',12,'FontName','TeX Gyre Schola Math'};
    font_axis14 = {'FontSize',14,'FontName','TeX Gyre Schola Math'};
    font_axis18 = {'FontSize',18,'FontName','TeX Gyre Schola Math'};
    font_axis22 = {'FontSize',22,'FontName','TeX Gyre Schola Math'};
    font_axis26 = {'FontSize',26,'FontName','TeX Gyre Schola Math'};
    font_latex18c = {'interpreter','latex','FontSize',18,'HorizontalAlignment','center','FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex8  = {'interpreter','latex','FontSize', 8,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex10 = {'interpreter','latex','FontSize',10,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex12 = {'interpreter','latex','FontSize',12,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex14 = {'interpreter','latex','FontSize',14,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex16 = {'interpreter','latex','FontSize',16,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex18 = {'interpreter','latex','FontSize',18,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex22 = {'interpreter','latex','FontSize',22,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex26 = {'interpreter','latex','FontSize',26,'FontName','TeX Gyre Schola Math','Color',[0 0 0]};
    font_latex = {'interpreter','latex'};
end

methods(Static)
    function ret = latexify_axis(ax, ax_fontsize)
        if isnumeric(ax)
            ax_fontsize = ax;
            ax = gca;
        end

        set(ax.Parent, 'Color', [1 1 1])

        set(ax, 'FontName','TeX Gyre Schola Math',...
            'GridColor', [0.1 0.1 0.1], 'MinorGridColor', [0.1 0.1 0.1],...
            'XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0]);

        set(ax,'TickLabelInterpreter', 'latex');

        xl = get(ax,'XLabel');
        xlFontSize = get(xl,'FontSize');
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', ax_fontsize)
        set(xl, 'FontSize', xlFontSize);

        yl = get(ax,'YLabel');
        ylFontSize = get(yl,'FontSize');
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', ax_fontsize)
        set(yl, 'FontSize', ylFontSize);

        zl = get(ax,'ZLabel');
        zlFontSize = get(zl,'FontSize');
        zAX = get(gca,'ZAxis');
        set(zAX,'FontSize', ax_fontsize)
        set(zl, 'FontSize', zlFontSize);

        if nargout > 0
            ret = ax;
        end
    end

    function ret = latexify_colorbar(cbar, fontsize)
        set(cbar, 'FontSize', fontsize, 'FontName','TeX Gyre Schola Math',...
            'TickLabelInterpreter', 'latex');
        cbar.Label.Color = [0 0 0];
        cbar.Label.Interpreter = 'latex';

        if nargout > 0
            ret = cbar;
        end
    end

    function ret_ = latexified_labels(ax, fontsize, xl, yl, zl)

        ret = cell(1,nargin - 2);
        if nargin > 2
            ret{1} = xlabel(ax,xl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);

        end

        if nargin > 3
            ret{2} = ylabel(ax,yl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);
        end

        if nargin > 4
            ret{3} = zlabel(ax,zl,'interpreter','latex','FontSize', fontsize,'FontName','TeX Gyre Schola Math','Color',[0 0 0]);
        end
        
        if nargout > 0
            ret_ = ret;
        end
    end

end

end
