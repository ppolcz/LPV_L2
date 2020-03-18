function [ret] = snp_publish(doch,event)
%% 
%  
%  file:   snp_publish.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.12.10. Saturday, 15:59:21
%
%%

G = pglobals;

% Get active editor
active = matlab.desktop.editor.getActive;
f = pcz_resolvePath(active.Filename);

% Load persistent object from the base workspace
try c = evalin('base','persist'); catch; c = []; end
persist = pcz_persist(f.path, c);
if isempty(c)
    assignin('base','persist',persist);
end

if ~exist(persist.pub_absdir,'dir')
    mkdir(persist.pub_absdir)
end

% Publish options
opt.format = 'html';
opt.outputDir = persist.pub_absdir;
opt.stylesheet = [G.ROOT '/publish.xsl' ];

% Publis + tidy html code
pub_output = publish(f.path, opt); 
system([ 'tidy -im ' pub_output])
fprintf('\ngenerated output: \n%s\n\n', pub_output)

% Copy html code to clipboard
persist.pub_html = sprintf('<a class="" href="<?php echo base_url(''%s/%s.html'') ?>">%s</a>', ...
    persist.pub_reldir, f.bname, f.bname);
clipboard('copy', persist.pub_html);

% Update persist object of the base workspace
persist.pub_output = pub_output;
assignin('base', 'persist', persist)

subl(pub_output)
open(pub_output)

end