% CellTracking (CT)
%
% Will run the C++ tracking from AIA interface.
%
% version 1.4
% Boris Guirao

%% Creates folder %%

if ~exist(trackingFolder,'dir')
    mkdir(trackingFolder);
end

        
%% Running C++ tracking %%

% Tracking patches:
CToptions = '1111 -T';
% NB: first set of 0-1 characters determines the patches used in the
% tracking; last “-T” determines whether the tracking will use
% transposed version of segmented and PIV images (also swapping u and v displacement matrices in that case).

% path to first transposed image:
startTransImagePath = [pathFolderRES filesep filename num2str(startFrame, digitsFormat) '.' imageFormat];
% path to first transposed PIV "u" file:
startPIVuImageFile = [pathFolderPIV filesep 'Results' filesep filenamePIV '_u_' num2str(startFrame,digitsFormat) '.png']; % 1.4
% startPIVuImageFile = [pathFolderPIV filesep PIVgrid 'Grid' filesep 'Results' filesep filenamePIV '_u_' num2str(startFrame,digitsFormat) '.png'];


cmd = [trackingExe ' ' startTransImagePath ' ' trackingFolder ' ' num2str(startFrame) ' ' num2str(finalFrame) ' ' startPIVuImageFile ' ' CToptions];
system(cmd);


%% History %%

% IMPROVEMENTS
% Nothing, it is perfects it is!

% 14/05/2018: 1.4
% - updated "startPIVuImageFile" now that "pathFolderPIV" contains the grid
% size subfolder.
% - removed "transFolder" and call to "pathFolderTRES" related to transpose
% of PIV images

% 07/05/2018: 1.3
% - stopped checking existence of last "correspondence" txt file here
% (moved it to SAP)

% 28/04/2018: 1.2
% - removed parts transposing segmented and PIV images as new tracking does NOT requires it anymore!
% - cleaned up commented parts 

% 30/03/2018: 1.1
% - change transpose segmentation image folder to fit new organisation

% - no "noPIV" mode right now (because NOT creating fake PIV png full of zeros!)

% **C++tracking** IMPROVEMENTS
% - parameter "noPIV" has to be passed as argument or option instead
% of writing png files with 0 everywhere that will be loaded like the
% regular ones!!

% 08/03/2018: creation

