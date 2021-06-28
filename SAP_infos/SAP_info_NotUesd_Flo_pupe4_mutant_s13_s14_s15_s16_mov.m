%% Frames to process %%

startFrame = 1;       % starts @ 1        
finalFrame = 182;       % ends @ 182

%% Path to raw images %%

pathFolderRaw = 'D:\Florencia\Analysis_Florencia_Big_RhoGEF4\Flo_pupe4_mutant_s13_s14_s15_s16_mov';
filenameRaw = {'Flo_pupe4_mutant_s13_s14_s15_s16_mov_'};                  % can be multiple: {'pten_myo_' ; 'pten_baz_' ; ...}. NB: PIV will use FIRST listed
signalName = {'Cad'};                       % "Cad", "Myo", "Mud"... (start with capital, names must be used consistently between animals)
nDigits = 4;                                % number of digits that WAS used for naming RAW IMAGES.
imageFormatRaw = 'tif';                     % image extension

%Cells below MUST BE changed/actualised


%% Synchronization & Scale %%

halfNotum = 'r';                            % Chose between 'l' (left) 'r' (right) and 'b' 'both'
frameRef = 49;                              % frame # corresponding to "timeRef" determined by "TimeRegistration"
timeRef = '20h15';                          % "timeRef" ('HH f:MM' or decimal) corresponding to "frameRef" above
dt = 5;                                     % time IN MINUTES between two frames
temperature = 25;                           % temperature (25 or 29C) at which development was filmed
yMid = [];                               % y of midline IN PIXELS (leave empty [] if unknown)
scale1D = 0.161;                            % Length of one pixel IN MICRONS (if 1 pixel is 0.1 um enter 0.1) 

%% Box/Grid specifications (PIV, CPT, AOS, TA, SM) %%

% NB: values obtained from rescaling of BIGwt2 alone on archetype including all eLife half-BIGs (20hAPF)
xFactor = 1;                                % X scaling factor from "RescaleAnimals". Animals LONGER than archetype: xFactor > 1
yFactor = 1;                                % Y scaling factor from "RescaleAnimals". Animals WIDER than archetype: yFactor > 1

% boxSize = [];  	% grid compartment size IN PIXELS [WIDTH HEIGHT]. One value => square grid
% xyStart = [];                               % xy coordinates (Upper Left Corner)of FIRST compartment. If empty, starts at image center.
% gridSize = [];                                      % number of compartments [ny nx]. If empty grid fills up image.


boxSize = [40 40]/scale1D  ./ [xFactor yFactor];  	% grid compartment size IN PIXELS [WIDTH HEIGHT]. One value => square grid
xyStart = [];                               % xy coordinates (Upper Left Corner)of FIRST compartment. If empty, starts at image center.
gridSize = [];                                      % number of compartments [ny nx]. If empty grid fills up image.

% boxSize = [128 128] ./ [xFactor yFactor];  	% grid compartment size IN PIXELS [WIDTH HEIGHT]. One value => square grid
% xyStart = [227 235];                       	% xy coordinates (Upper Left Corner)of FIRST compartment. If empty, starts at image center.
% % xyStart = [225 225];                     	% xy coordinates (Upper Left Corner)of FIRST compartment. If empty, starts at image center.
% gridSize = [28 28];                           % number of compartments [ny nx]. If empty grid fills up image.
% % gridSize = [14 14];                         % number of compartments [ny nx]. If empty grid fills up image.

% JESUS
% boxSize = [60 128]./[xFactor yFactor]/scale1D;     % 60 um x 128 um box in pixels
% dist2Macro1 = 65/scale1D;                           % 65 um in pixels
% xMacro1 = 324;                                       % 1st macro X coord clicked @24 hAPF
% xyStart = [xMacro1+dist2Macro1 yMid-boxSize(2)/2];  % defining TOP LEFT corner of box. NB: yMid IS REQUIRED
% gridSize = [1 1];
        

%% Path to clone Mask %%

%cloneMaskFrame = 9;                    % frame number corresponding to clone mask
% cloneName = 'COMP';                         % name to give to clone ("a-act", "fat"...)
% path2cloneMask = 'D:\BigMovies\BIGwt2\Masks\COMP_mask_BIGwt2_0009.png'; 
% cloneName = 'DEL';                         % name to give to clone ("a-act", "fat"...)
% path2cloneMask = 'D:\BigMovies\BIGwt2\Masks\DEL_mask_BIGwt2_0009.tif'; 

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

