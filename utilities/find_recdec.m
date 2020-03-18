function [ret1,ret2,ret3] = find_recdec(theFraction, varargin)
%% find_recdec
%  
%  File: find_recdec.m
%  Directory: utilities
%  
%  Created on 2019. June 07. (2019a)
%
%  Sources: https://www.mathworks.com/matlabcentral/answers/249333-code-to-recognize-repeating-numbers
%  Answer by Image Analyst (https://www.mathworks.com/matlabcentral/profile/authors/1343420-image-analyst)
% 
%  Improved by Peter Polcz (ppolcz@gmail.com) 
% 

%%
 
args.tolerance = 1e-10;
args.symbolic = 1;
args = parsepropval(args,varargin{:});

%%

min_denominator = 2;
max_denominator = 1000;
step = 1000;
max_nr_of_iterations = 2;

ret1 = theFraction;
ret2 = 1;

integer_part = floor(theFraction);

theFraction = theFraction - integer_part;

% Handle integer numbers
if theFraction < args.tolerance
    ret1 = integer_part;
    if args.symbolic == 1
        ret2 = 1;
        ret3 = theFraction;
    else
        ret2 = theFraction;
    end
    return
end

for i = 1:max_nr_of_iterations

    denominators = min_denominator:max_denominator;
    
    % Find potential numerators, rounded to the closest integer:
    numerators = round(theFraction * denominators);
    
    % We need to get rid of zeros.
    zeroIndexes = numerators == 0;
    numerators(zeroIndexes) = [];
    denominators(zeroIndexes) = [];
    
    % Now get the ratios of integers.
    ratios = numerators ./ denominators;
    differences = abs(ratios - theFraction);

    % Find the min difference:
    [minDifference, indexOfMin] = min(differences);

    if minDifference > args.tolerance
        min_denominator = max_denominator;
        max_denominator = max_denominator*step;
        continue
    end
    
    if args.symbolic == 1
        ret1 = sym(numerators(indexOfMin)) / sym(denominators(indexOfMin)) + sym(integer_part);
        ret2 = minDifference;
    else
        ret1 = integer_part*denominators(indexOfMin) + numerators(indexOfMin);
        ret2 = denominators(indexOfMin);
        ret3 = minDifference;
    end
    
    break
end


end


function standalone
%%

%%

theFraction = .31454545;  % Whatever number you're examining


format long g;
format compact;
fontSize = 20;

% If the fraction is within this, say it's a match
tolerance = 0.1;

denominators = 2:1000;
% Find potential numerators, rounded to the closest integer:
numerators = round(theFraction * denominators);
% We need to get rid of zeros.
zeroIndexes = numerators == 0;
numerators(zeroIndexes) = [];
denominators(zeroIndexes) = [];
% Now get the ratios of integers.
ratios = numerators ./ denominators;
differences = abs(ratios - theFraction);
subplot(2,1,1);
plot(ratios);
title('Ratios', 'FontSize', fontSize);
grid on;
subplot(2,1,2);
plot(differences);
title('Differences', 'FontSize', fontSize);
grid on;
% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 
% Find the min difference:
[minDifference, indexOfMin] = min(differences)
% Print the fraction that is closest:
message = sprintf('%.8f is approximately equal to %d / %d, which equals %.9f.\nThe difference is %.9f.',...
theFraction, round(numerators(indexOfMin)), denominators(indexOfMin), ...
numerators(indexOfMin)/denominators(indexOfMin), minDifference)
uiwait(helpdlg(message));

end