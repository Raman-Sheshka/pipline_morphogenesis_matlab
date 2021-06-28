%% Frames to process %%

startFrame = 9;         % starts @ 9             
finalFrame = 240;       % ends @ 240

%% Path to raw images %%

pathFolderRaw = 'C:\Users\ines\Documents\Python_Related\Evaluation_Stardist\BIGs\BIGwt2r';   
filenameRaw = {'BIGwt2r_'};      
signalName = {'Cad'};                       % "Cad", "Myo", "Mud"... (USE UPPER CASE)
nDigits = 4;                                % number of digits that WAS used for naming RAW IMAGES.
imageFormatRaw = 'tif';                     % image extension

%% Synchronization & Scale %%

halfNotum = 'r';                            % Chose between 'l' (left) 'r' (right) and 'b' 'both'
frameRef = 96;                              % frame # corresponding to "timeRef" determined by "TimeRegistration"
timeRef = '20h15';                          % "timeRef" ('HH:MM' or decimal) corresponding to "frameRef" above
dt = 5;                                     % time IN MINUTES between two frames
temperature = 25;                           % temperature (25 or 29C) at which development was filmed
scale1D = 0.322;                            % Length of one pixel IN MICRONS (if 1 pixel is 0.1 um enter 0.1) 

%% Box/Grid specifications (PIV, CPT, AOS, TA, SM) %%

% JESUS (rescaling at 26 hAPF)
% yMid = 76;                                 % y of midline IN PIXELS (leave empty [] if unknown)
% xFactor = 1.00;                                % X scaling factor from "RescaleAnimals". Animals LONGER than archetype: xFactor > 1
% yFactor = 0.97;                                % Y scaling factor from "RescaleAnimals". Animals WIDER than archetype: yFactor > 1
% xMacro1 = 318;                                       % 1st macro X coord clicked @26 hAPF
% boxSize = [60 yMid*scale1D+64]./[xFactor yFactor]/scale1D;     % 60 um x 128 um box in pixels
% dist2Macro1 = 65/scale1D;                           % 65 um in pixels
% xyStart = [xMacro1+dist2Macro1 1];  % defining TOP LEFT corner of box. NB: yMid IS REQUIRED
% gridSize = [1 1];

% rescaling at 20 hAPF
yMid = 84;                                 % y of midline IN PIXELS (leave empty [] if unknown)
xFactor = 0.849;                                % X scaling factor from "RescaleAnimals". Animals LONGER than archetype: xFactor > 1
yFactor = 0.829;                                % Y scaling factor from "RescaleAnimals". Animals WIDER than archetype: yFactor > 1
boxSize = [40 40]/scale1D  ./ [xFactor yFactor];   % grid compartment size IN  IN MICRON. One value => square grid
xyStart = [715 417];                        % xy coordinates (Upper Left Corner)of FIRST compartment. If empty, starts at image center.
gridSize = [];                              % number of compartments [ny nx]. If empty grid fills up image.
        
% AURELIEN
% Animal = 'BIGwt2r';
% path2CTDmacroFile = ['D:\BigMovies\' Animal '\SAP_' Animal '\CTD_' Animal '\macroCells.mat'];
% path2SIAbu = ['D:\BigMovies\' Animal '\SAP_' Animal '\SIA_' Animal '\Backups'];
% MakeAurelienGrid;
        

%% Path to clone Mask %%

cloneName = '';                         % name to give to clone ("a-act", "fat"...)
cloneMaskFrame = [];                    % frame number corresponding to clone mask
path2cloneMask = ''; 

%% Scale Bar & Display %%

scaleBarLength = 50;                    % scale bar length in microns
scaleBarWidth = 5;                      % scale bar width (in PIXELS). Used 3 for BIG
xyOffset = [30 30];                     % scale bar offset @ bottom right of image IN PIXELS ([30 30] for Potts)
colorBarXYWH = [0.75 0.025 0.15 0.02];    % [left, bottom, width, height] in figure units
imageFading = 0.7;                      % fading of background image for tensor maps (0.7 for eLife article)
fontSizeInfo = 25;                      % 27 TRBL8(res 200), 35 BIG_X,TRBL9,TRBL4(res 150), 55 TRBL7(res 100), 20 BIG (res 300), 40 JESUS(res 200), but 27 (res300) for cad4
imageResolution = 200;                  % 300 used for formalism article, 200 for Potts CppTD images, between 300 and 150 for CppTD images

%% "SAP_parameters" call %%

Animal = mfilename;
Animal = Animal(10:end);
SAP_parameters

