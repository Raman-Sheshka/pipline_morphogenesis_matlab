%% Frames to process %%

startFrame = 1;         % starts @ 1           
finalFrame = 216;       % ends @ 216 
%% Path to raw images %%

pathFolderRaw = 'E:\Melanie_Analyses_in_progress\AnalysesDivisions\20210131_cad3xGFP_ctle_pupe2';
filenameRaw = {'Projected_notum_ctle_cad3xGFP_40x_bin1_1_w2CSU-488_'};                  % can be multiple: {'pten_myo_' ; 'pten_baz_' ; ...}. NB: PIV will use FIRST listed
signalName = {'Cad'};                       % "Cad", "Myo", "Mud"... (start with capital, names must be used consistently between animals)
nDigits = 4;                                % number of digits that WAS used for naming RAW IMAGES.
imageFormatRaw = 'tif';                     % image extension

%% Synchronization & Scale %%

halfNotum = 'b';                            % Chose between 'l' (left) 'r' (right) and 'b' 'both'
frameRef = 1;                              % frame # corresponding to "timeRef" determined by "TimeRegistration"
timeRef = '16h30';                          % "timeRef" ('HH:MM' or decimal) corresponding to "frameRef" above
dt = 5;                                     % time IN MINUTES between two frames
temperature = 29;                           % temperature (25 or 29C) at which development was filmed
yMid = [935];                                 % y of midline IN PIXELS (leave empty [] if unknown)
scale1D = 0.161;                            % Length of one pixel IN MICRONS (if 1 pixel is 0.1 um enter 0.1)                         % Length of one pixel IN MICRONS (if 1 pixel is 0.1 um enter 0.1) 

%% Box/Grid specifications (PIV, CPT, AOS, TA, SM) %%

xFactor = 1;                                % X scaling factor from "RescaleAnimals". Animals LONGER than archetype: xFactor > 1
yFactor = 1;                                % Y scaling factor from "RescaleAnimals". Animals WIDER than archetype: yFactor > 1
boxSize = [40 40]/ scale1D ./ [xFactor yFactor];   % grid compartment size IN PIXELS [WIDTH HEIGHT]. One value => square grid
xyStart = [1 1];                        % xy coordinates (Upper Left Corner)of FIRST compartment. If empty, starts at image center.
gridSize = [];                              % number of compartments [ny nx]. If empty grid fills up image.
        

%% Path to clone Mask %%

% cloneName = 'DIAP1_1a';                         % name to give to clone ("a-act", "fat"...)
% cloneMaskFrame = [1];                    % frame number corresponding to clone mask
% path2cloneMask = 'D:\Maria\DIAP_c1_mov1a\Clone_Mask\LastMask.tif'; 

%% Scale Bar & Display %%

scaleBarLength = 50;                    % scale bar length in microns
scaleBarWidth = 3;                      % scale bar width (in PIXELS). Used 3 for BIG
xyOffset = [30 70];                     % scale bar offset @ bottom right of image IN PIXELS ([30 30] for Potts)
colorBarXYWH = [0.75 0.025 0.15 0.02];    % [left, bottom, width, height] in figure units
imageFading = 0.7;                      % fading of background image for tensor maps (0.7 for eLife article)
fontSizeInfo = 18;                      % 27 TRBL8(res 200), 35 BIG_X,TRBL9,TRBL4(res 150), 55 TRBL7(res 100), 20 BIG (res 300), 40 JESUS(res 200), but 27 (res300) for cad4
imageResolution = 200;                  % 300 used for formalism article, 200 for Potts CppTD images, between 300 and 150 for CppTD images

%% "SAP_parameters" call %%

Animal = mfilename;
Animal = Animal(10:end);
SAP_parameters

