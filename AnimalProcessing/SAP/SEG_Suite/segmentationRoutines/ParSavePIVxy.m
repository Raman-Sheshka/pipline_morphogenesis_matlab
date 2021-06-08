function ParSavePIVxy(path, filename, x, y)  %#ok<INUSD>

% Version 1.0
% Boris Guirao

% Workaround to be able to save the 2 variables x and y containing the
% locations where displacement fields u v are calculated.

% NB: THIS FUNCTION IS VARIABLE SPECIFIC AND ONLY WORKS FOR "x" and "x"

full_path = fullfile(path, [filename '.mat']);
save(full_path, 'x', 'y')

%% History

% 17/06/2010: creation

