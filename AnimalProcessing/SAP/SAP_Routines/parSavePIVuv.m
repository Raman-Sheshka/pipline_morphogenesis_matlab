function parSavePIVuv(path, filename, digitsFormat, index, u, v, uq, vq)  %#ok<INUSD>
%
% parSavePIVuv(path, filename, digitsFormat, index, u, v, uq, vq)
%
% Workaround to be able to save the 4 variables u,v and uq,vq containing the displacement fields at each iteration of the parfor loop for velocity
% fields to be used with segmentation/tracking (u,v) and for quantitative analysis (uq,vq).
%
% NB: THIS FUNCTION IS VARIABLE SPECIFIC AND ONLY WORKS FOR "u","v","uq", "vq"
%
% Version 1.2
% Boris Guirao

%% Code %%

fullPath = fullfile(path, [filename num2str(index, digitsFormat) '.mat']);
save(fullPath, 'u', 'v', 'uq', 'vq')

%% History

% 26/02/2014: 1.2
% - now saves uq,vq, hence only compatible with "Particle_Image_Velocimetry" v2.2+

% 21/06/2010: 1.1
% - now saves u,v

% 17/06/2010: creation

