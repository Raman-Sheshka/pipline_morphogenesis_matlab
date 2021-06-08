% CreateArchetype
%
% This program gathers the positions of different landmarks clicked in different movies and uses them to
% create an archetype animal
%
% NB: Rule is Posterior to the left side of the image, anterior to the right, left at the top of image
%                 Left
%                   |
%                   |
%                   |
% Posterior ---------------- Anterior
%                   |
%                   |
%                   |
%                 Right
%
% Rigaud Stephane
% Boris Guirao
% Anais Bailles
% Isabelle Bonnet

clear; close all; clc
version = '2';



%% PARAMETERS %%

% inputFolder = uigetdir();
inputFolder = 'D:\JESUS\Actn\mf3\STR_mf3\spaceReg_mf3_47';
% inputFolder = 'D:\JESUS\mfX_rescale';
timeRef = '1840';
nMacroMin = 3;                % MIN number of macrochaetes AVAILABLE (ON EACH SIDE OF ANIMAL if halfNotum = 'b') in all animals to process (USED TO CALCULATE BARYCENTERS)
halfNotum = 'b';
scale1D = 0.216;

% inputFolder = 'D:\BigMovies\spatial_registration\archetypeElife_26';
% timeRef = '2600';
% inputFolder = 'D:\archetype_full\elife'; % Folder containing all the txt files of Macro and Midline coordinates to be used for space registration
% inputFolder = 'D:\archetype_full\maria'; % Folder containing all the txt files of Macro and Midline coordinates to be used for space registration
% nMacroMin = 5;                % MIN number of macrochaetes AVAILABLE (ON EACH SIDE OF ANIMAL if halfNotum = 'b') in all animals to process (USED TO CALCULATE BARYCENTERS)
% NB: if nMacroMin = 5 and halfNotum = 'b', 10 macro will actually be used, thereby respecting the symmetry of the animals
% halfNotum = 'r';
% scale1D = 0.322;

originType = 2;               % 2 to use the barycentre of the macrochaetae as the coordinate origin, 1 to use the y midline, 0 to use the first macrochaetae
makeSVG = false;


nLandmarksMax = 8;
if halfNotum == 'b'
    nLandmarksMax = 16;
end



%% Loading files from inputFolder
[moviesName,dataMacroX,dataMacroY,dataMidLin,dataNeck] = LandmarkRescaleLoader(inputFolder, halfNotum);
nMovies = size(moviesName,1);




%% Switching scale, from image to archetype scale

% determine if mid landmarks play in the origins
flagMidlandmarks = false;
%--------------------------------------------------------------------------
% % Use of the mid point into the origin calculation
% if all( ~isnan(dataNeck(3,:)) & ~isnan(dataMidLin(3,:)))
%   flagMidlandmarks = true;  
%   fprintf('warning: Mid landmarks will be use for origin calculation.\n');
% end
%--------------------------------------------------------------------------

[originX, originY] = ComputeOrigin(originType, dataMacroX, dataMacroY, dataMidLin, dataNeck, ...
                                   nMacroMin, halfNotum, flagMidlandmarks);

% Coordinates refered to origin
oDataMacroX = dataMacroX - repmat(originX,[nLandmarksMax 1]); % (pxl)
oDataMacroY = dataMacroY - repmat(originY,[nLandmarksMax 1]); % (pxl)
oDataMidlin = dataMidLin - repmat(originY,[3 1]);             % (pxl)
oDataNeck   = dataNeck   - repmat(originX,[3 1]);             % (pxl)






%% Rescaling

outputFolder = [inputFolder filesep 'archetype_' timeRef '_nMacroUsed=' num2str(nMacroMin)]; % 1.0
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end

nLandmarksMax = nLandmarksMax + 1;

% Pixels scales matrix
scaleMat = repmat(scale1D, nLandmarksMax, nMovies);

Xall = [oDataMacroX ; oDataNeck(3,:)];   % add Neck coordinate (can be NaN)
Xall = Xall .* scaleMat;            % Convert pixel coordinate into um

Yall = [oDataMacroY ; oDataMidlin(3,:)]; % add Midline coordinate (can be NaN)
Yall = Yall .* scaleMat;            % Convert pixel coordinate into um


[moviesA, moviesB, scalingFactorX, scalingFactorY] = ComputeRescale(Xall, Yall, halfNotum);


% Standard deviation around archetype position of Archetype and Movies to be rescale;
stdMovies = [nanstd(Xall,0,2) nanstd(Yall,0,2)];

% Standard deviation before and after rescaling; ** Use of Xall and Yall (1.3) **
std_not_rescaled = [ nanstd(Xall, 0, 2)  nanstd(Yall, 0, 2) ];
std_rescaled     = [ nanstd(Xall .* repmat(scalingFactorX,nLandmarksMax,1),0,2) ...
    nanstd(Yall .* repmat(scalingFactorY,nLandmarksMax,1),0,2) ];



%% Save outputs %%

macroRawXY = [ dataMacroX ; dataMacroY ];

% rescaling output information
rescalingOutput = cell(8,nMovies+1);
rescalingOutput(1,1)       = {'Name'};
rescalingOutput(1,2:end)   = moviesName;
rescalingOutput(2:3,1)     = [{'O_x'}; {'O_y'}];
rescalingOutput(2:3,2:end) = num2cell([originX ; originY]);
rescalingOutput(4:5,1)     = [{'factor_x'}; {'factor_y'}];
rescalingOutput(4:5,2:end) = num2cell([scalingFactorX ; scalingFactorY]);
rescalingOutput(6:7,1)     = [{'Y_ML'}; {'X_NK'}];
rescalingOutput(6:7,2:end) = num2cell([dataMidLin(3,:) ; dataNeck(3,:)]);
rescalingOutput(8,1)       = {'scale_1D'};
rescalingOutput(8,2:end)   = num2cell(scaleMat(1,:));
% write output in .mat and .txt format
save([outputFolder filesep 'rescalingOutput_' num2str(originType)],'rescalingOutput');
dlmcell([outputFolder filesep 'rescalingOutput_' num2str(originType) '.txt'],rescalingOutput,' ')

% std ouput information
stdOutput = cell(6,2);
stdOutput(1:2,1) = [{'std_Movies_x'} ; {'std_Movies_y'}];
stdOutput(1:2,2) = num2cell([nanmean(stdMovies(:,1)) ; nanmean(stdMovies(:,2))]);
stdOutput(3:4,1) = [{'std_not_rescaled_x'} ; {'std_not_rescaled_y'}];
stdOutput(3:4,2) = num2cell([nanmean(std_not_rescaled(:,1)) ; nanmean(std_not_rescaled(:,2))]);
stdOutput(5:6,1) = [{'std_rescaled_x'} ; {'std_rescaled_x'}];
stdOutput(5:6,2) = num2cell([nanmean(std_rescaled(:,1)) ; nanmean(std_rescaled(:,2))]);
% write output in .mat and .txt format
save([outputFolder filesep 'stdOutput_' num2str(originType)],'stdOutput');
dlmcell([outputFolder filesep 'stdOutput_' num2str(originType) '.txt'],stdOutput,' ')




%% Graphical Outputs %%

cmap = colormap(jet); close;
factor = size(cmap,1) / nMovies - 1;
size1 = 1:nLandmarksMax-1;
size2 = [1:nLandmarksMax-1 1:nLandmarksMax-1];

% determination of max values for plot (1.2,1.3)
Xmin = min(Xall(:));
Xmax = max(Xall(:));
Ymax = max(abs(Yall(:)));
Ymin = -Ymax; % for symmetry % horizontal axis
if ~strcmp(halfNotum,'b')
    Ymin = min(Yall(:));
end
% introduces some leeway with axis box
Xmax = 1.2 * Xmax; % 1.2
Xmin = 1.2 * Xmin; % 1.2
Ymin = 1.1 * Ymin; % 1.1
Ymax = 1.1 * Ymax; % 1.1

% Raw movies vs archetype
figure(1)
hold on
axis equal                  % 1.2
axis ij                     % 1.2
box on                      % 1.2
axis([Xmin Xmax Ymin Ymax]) % 1.2
for n = 1:nMovies
    % color code
    color = cmap(round(factor * n),:);
    % scaling
    midlineVectorY = [oDataMidlin(n) oDataMidlin(n)] .* scaleMat(1,n); % (um)
    neckVectorX    = [oDataNeck(n) oDataNeck(n)]     .* scaleMat(1,n); % (um)
    macroVectorX   = oDataMacroX(:,n)                .* scaleMat(1,n); % (um)
    macroVectorY   = oDataMacroY(:,n)                .* scaleMat(1,n); % (um)
    % midline plot
    line([Xmin Xmax], midlineVectorY, 'LineStyle', ':', 'Color', color, 'LineWidth', 2)
    % macrochaete plot
    scatter(macroVectorX, macroVectorY, 20, color, 'LineWidth', 2);
    % neck plot
    scatter(neckVectorX, midlineVectorY, 20, color, 'LineWidth', 2);
end
% archetype plot
plot(moviesA(1:end-1,1), moviesA(1:end-1,2), 'xb', 'LineWidth', 2, 'MarkerSize', 10) % machochaete
plot(moviesA(end,1), moviesA(end,2), 'xb', 'LineWidth', 2, 'MarkerSize', 5)          % neck
line([Xmin Xmax], [moviesA(end,2) moviesA(end,2)], 'Color', [0 0 0], 'LineWidth', 2) % midline
% origin (0,0)
plot(0, 0, '*k', 'LineWidth', 2, 'MarkerSize', 5)
% save figure
title('Raw Movies vs Archetype')
print('-dpng', '-r400', [outputFolder filesep 'raw_movies_vs_archetype_' num2str(originType) '.png']); % 1.1
if makeSVG
    plot2svg([outputFolder filesep 'raw_movies_vs_archetype_' num2str(originType) '.svg'], figure(1), 'png'); % 1.1
end
close


% Rescaled movies vs archetype
figure(2)
hold on
axis equal                  % 1.2
axis ij                     % 1.2
box on                      % 1.2
axis([Xmin Xmax Ymin Ymax]) % 1.2
for n = 1:nMovies
    % color code
    color = cmap(round(factor * n),:);
    % scaling and registering
    midlineVectorY = [oDataMidlin(n) oDataMidlin(n)] .* scalingFactorY(n) .* scaleMat(1,n); % (um)
    neckVectorX    = [oDataNeck(n)   oDataNeck(n)]   .* scalingFactorY(n) .* scaleMat(1,n); % (um)
    macroVectorX = scalingFactorX(n) .* oDataMacroX(:,n) .* scaleMat(1,n); % (um)
    macroVectorY = scalingFactorY(n) .* oDataMacroY(:,n) .* scaleMat(1,n); % (um)
    % mideline plot
    line([Xmin Xmax], midlineVectorY, 'LineStyle', ':', 'Color', color, 'LineWidth', 2)
    % macrochaete plot
    scatter(macroVectorX, macroVectorY, 20, color, 'LineWidth', 2)
    % neck plot
    scatter(neckVectorX, midlineVectorY, 20, color, 'LineWidth', 2);
end
% archetype plot
plot(moviesA(1:end-1,1), moviesA(1:end-1,2), 'xb', 'LineWidth', 2, 'MarkerSize', 10) % machochaete
plot(moviesA(end,1), moviesA(end,2), 'xb', 'LineWidth', 2, 'MarkerSize', 5)          % neck
line([Xmin Xmax], [moviesA(end,2) moviesA(end,2)], 'Color', [0 0 0], 'LineWidth', 2) % midline
% origin (0,0)
plot(0, 0, '*k', 'LineWidth', 2, 'MarkerSize', 5)
% save figure
title('Rescale Movies vs Archetype')
print('-dpng', '-r400', [outputFolder filesep 'rescale_movies_vs_archetype_' num2str(originType) '.png']); % 1.1
if makeSVG
    plot2svg([outputFolder filesep 'rescale_movies_vs_archetype_' num2str(originType) '.svg'], figure(1), 'png'); % 1.1
end
close




%%  Saving txt file indicating date and version used in "Save_folder" (1.3) %%

today = datestr(now,29);                      % format 29 displays date yyyy-mm-dd style. Look up date for details
txtFilename = [today '_Archetype_' version '.txt'];
% Writing main parameters in txt file (3.1)
parameterCell = {   'Main Parameters:',[];
    [],[];
    'inputFolder = ', inputFolder;
    'nMacroMin = ', nMacroMin;
    'nMovies = ', nMovies;
    'scale1D = ', scale1D;
    'halfnotum = ', halfNotum;
    'originType = ', originType;
    'useMidPoint = ', flagMidlandmarks};
dlmcell([outputFolder filesep txtFilename], parameterCell,' ');




%% Generate final archetype right and full %%

[Full, Half] = CreateSymmetricArchetype(moviesA(1:end-1,1),moviesA(1:end-1,2), ... 
                                        moviesA(end,:), halfNotum, nMacroMin, flagMidlandmarks);
Full.scale1D = scale1D;
Half.scale1D = scale1D;

originShift = Half.origin;
if halfNotum == 'r'
originShift = Full.origin;
end

Xmin = min(min(Full.macrochaete(:,1)),Xmin) * 1.2;
Xmax = max(max(Full.macrochaete(:,1)),Xmax) * 1.2;
Ymin = min(min(Full.macrochaete(:,2)),Ymin) * 1.1;
Ymax = max(max(Full.macrochaete(:,2)),Ymax) * 1.1;

figure(3)
hold on
axis equal                  % 1.2
axis ij                     % 1.2
box on                      % 1.2
axis([Xmin Xmax Ymin Ymax]) % 1.2
for n = 1:nMovies
    % we apply the shifted origin
    midlineVectorY = oDataMidlin(n)   .* scaleMat(1,n) - Full.origin(2); % (um)
    neckVectorX    = oDataNeck(n)     .* scaleMat(1,n) - Full.origin(1); % (um)
    macroVectorX   = oDataMacroX(:,n) .* scaleMat(1,n) - Full.origin(1); % (um)
    macroVectorY   = oDataMacroY(:,n) .* scaleMat(1,n) - Full.origin(2); % (um)
    scatter(macroVectorX, macroVectorY, 20, color, 'LineWidth', 2);
    scatter(neckVectorX, midlineVectorY, 20, color, 'LineWidth', 2);
end
% archetype plot
plot(Full.macrochaete(:,1), Full.macrochaete(:,2), 'xb', 'LineWidth', 2, 'MarkerSize', 10) % macrochaete
plot(Full.midpoint(1), Full.midpoint(2), 'xb', 'LineWidth', 2, 'MarkerSize', 10)           % neck
line([Xmin Xmax], [Full.midpoint(2) Full.midpoint(2)], 'Color', [.8 .8 .8]);               % midline
% origin (0,0)
plot(0, 0, '*k', 'LineWidth', 2, 'MarkerSize', 5)
hold off
% save figure
title('Full symetric archetype')
print('-dpng', '-r400', [outputFolder filesep 'full_symetric_archetype_' num2str(originType) '.png']); % 1.1
if makeSVG
    plot2svg([outputFolder filesep 'full_symetric_archetype_' num2str(originType) '.svg'], figure(3), 'png'); % 1.1
end
close

figure(3)
hold on
axis equal                  % 1.2
axis ij                     % 1.2
box on                      % 1.2
axis([Xmin Xmax Ymin Ymax]) % 1.2
for n = 1:nMovies
    % we apply the shifted origin
    midlineVectorY = oDataMidlin(n)   .* scaleMat(1,n) - Half.origin(2); % (um)
    neckVectorX    = oDataNeck(n)     .* scaleMat(1,n) - Half.origin(1); % (um)
    macroVectorX   = oDataMacroX(:,n) .* scaleMat(1,n) - Half.origin(1); % (um)
    macroVectorY   = oDataMacroY(:,n) .* scaleMat(1,n) - Half.origin(2); % (um)
    scatter(macroVectorX, macroVectorY, 20, color, 'LineWidth', 2);
    scatter(neckVectorX, midlineVectorY, 20, color, 'LineWidth', 2);
end
% archetype plot
plot(Half.macrochaete(:,1), Half.macrochaete(:,2), 'xb', 'LineWidth', 2, 'MarkerSize', 10)
plot(Half.midpoint(1),Half.midpoint(2), 'xb', 'LineWidth', 2, 'MarkerSize', 10)
line([Xmin Xmax], [Half.midpoint(2) Half.midpoint(2)],'Color',[.8 .8 .8]);
% origin (0,0)
plot(0, 0, '*k', 'LineWidth', 2, 'MarkerSize', 5)
% save figure
hold off
title('Half symetric archetype')
print('-dpng', '-r400', [outputFolder filesep 'half_symetric_archetype_' num2str(originType) '.png']); % 1.1
if makeSVG
    plot2svg([outputFolder filesep 'half_symetric_archetype_' num2str(originType) '.svg'], figure(3), 'png'); % 1.1
end
close

% save archetype
save([outputFolder filesep 'FullArchetype'],'-struct','Full');
save([outputFolder filesep 'HalfArchetype'],'-struct','Half');




%% History %%

% IMPROVEMENTS:
% - support cases where initial clicks starts with NaNs
% - support pupae oriented along vertical axis
% - load scale1D of each animal by loading its AIA_info file (currently only one scale1D authorized for all movies)
% - add symmetrie evaluation during archetype calculation

% 13/10/2016: 2 (Stephane)
% - PIPELINE MODIFICATION: Now rescale is performed in two separated steps. First step is the creation of an archetype as a set of similar animal
% rescaled and mean together to form one archetype animal using the process CreateArchetype. Second step is, using RescaleAnimals, the rescale of N animal on a selected
% already generated archetype.
% - two input folder: movie folder and archetype folder
% - remove of the archetype calculation. The process now only rescale movies on an archetype already calculated
% - add the mid point defined by the Y of the midline and the X of the neck. Used in the calculation of the rescaling factors.
% - (to be removed) add the mid point into the origin calculation.

% 20/07/2016: 1.3 (Boris)
% - FIXED MISTAKE: previous code was NOT using "okRowsWT" to only limit calculations to available macro in WT archetype, sometime leading to
% wrong scalingFactorX/Y!!
%
% NB1: "okRowsWT" also indicates rows of WT_A with macro coordinates: THIS IS THE MAXIMAL SET OF MACROCAETAE TO CONSIDER IN THE RESCALING (SEE
% eLIFE EQ. 10). IF an animal (WT or Mutant) has MORE macro available, the numerator of eq 10 will NOT take them into account (NaNs in WT_A
% giving NaN in the product) BUT THE DENOMINATOR WILL, thereby yielding an erroneous value of alpha for this animal.
% Inversely, if an animal has LESS macro available, the non-NaN values in WT_A will be overriden by the NaN in the animal, thereby
% considering the same macro at numerator and denominator of eq 10.
%
% NB2: this bug was noticed when rescaling diap1UP1scut with wt2NEW
% - saving a txt file indicating the version of the program used and the parameters values

% 15/07/2016: 1.2 (Boris)
% - much better support of animals with both sides: "nMacroMin" will be used ON BOTH SIDES of the animal for the rescaling
% - now checks that the nMacroMin entered by user is actually met in all txt files (namely that there is no NaN listed)
% - automatically detects "nLandmarksMax" that is no longer a parameter
% - now properly extracts animal names from "Macro_coordinates_animal_side_b.txt" files
% - updated "ClickLandmars" so that 16 Macrochaetes should be clicked (instead of 14 before)
% - improved plots
% - thoroughly used "filesep" for mac compatibility

% 27/05/2016: 1.1 (Boris)
% - fixed mismatch between saved plot names and what was actually plotted

% 02/05/2016: Boris: changed name to "RescaleAnimals" (1.0) from "Scatter_landmarks_and_rescale"
% - renamed many variables
% - output folder now is a subfolder of input folder and is named after the number of macrochaetes used for registration (_nMacroUsed=)
% - new parameter "makeSVG" to avoid always saving them


% 09/2014: Ana?s
% /!\ the macrochaetae have to be picked up in a definite order

% 05/2012: Isabelle
