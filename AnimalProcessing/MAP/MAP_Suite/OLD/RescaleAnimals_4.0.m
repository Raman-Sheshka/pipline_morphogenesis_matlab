% RescaleAnimals
%
% * if rescaleMode = "Archetype":
% This program gathers the positions of different landmarks clicked in
% different animals and uses them to CREATE an archetype.
% NB: replaces "CreateArchetype" as of 4.0
%
% * if rescaleMode = "Rescale":
% This program gathers the positions of different landmarks clicked in
% different animals and RESCALES them on an archetype generated by
% "CreateArchetype".
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
% Stephane Rigaud
% Boris Guirao
% Anais Bailles
% Isabelle Bonnet
version = '4.0';


%% Loading files from each animal SR folder (3.0)%%

% Matrix of animals' scale1D (3.0)
scale1Dvec = NaN(nAvgAnimals,1);
dtVec = NaN(nAvgAnimals,1);                 % NOT USED SO FAR
temperatureVec = NaN(nAvgAnimals,1);        % NOT USED SO FAR
halfNotumVec = repmat('0',nAvgAnimals,1);   % creates char vector

% Defining temp quantities to avoid overwritting
AnimalAOA = Animal;
PLOTAOA = PLOT;
temperatureAOA = temperature;

for a = 1:nAvgAnimals

    % Loading this animal info (3.0)
    aName = avgAnimals{a};
    AIA_info_file = [AIAFolderName filesep 'AIA_info_' aName '.m' ];
    if exist(AIA_info_file,'file')
        run(AIA_info_file);
    else
        disp(['WARNING: "AIA_info_' aName '.m" not found in specified location! Trying Matlab saved paths...'])
        eval(['AIA_info_' aName]);
    end
    
    % Filling vectors
    scale1Dvec(a) = scale1D;
    dtVec(a) = dt;                      % NOT USED SO FAR
    temperatureVec(a) = temperature;    % NOT USED SO FAR
    halfNotumVec(a) = halfNotum;
    
    if a == 1
        dataMacroX = NaN(nMacroMAXspace,nAvgAnimals);
        dataMacroY = NaN(nMacroMAXspace,nAvgAnimals);
        dataMidLin = NaN(3,nAvgAnimals);
        dataNeck   = NaN(3,nAvgAnimals);
    end
    
    %%% Loading click txt files
    % Defining this animal "clickFrame" corresponding to "clickTimeAll"
    clickFrame = round(time2frame(clickTimeAll, timeRef, frameRef, dt));
    pathFolderSR =  [pathFolderSTR filesep 'spaceReg_'  Animal '_' num2str(clickFrame)];    % "SR" for "SpaceRegistration" (7.12)
    % Loading txt files for this animal
    [~,dataMacroX(:,a),dataMacroY(:,a),dataMidLin(:,a),dataNeck(:,a)] = LandmarkRescaleLoader(pathFolderSR, halfNotum);
end

% Checks all animals' "halfNotum" are the same
halfNotum = halfNotumVec(1);
allHalfNotumIdenticalTF = all(ismember(halfNotumVec,halfNotum),1);
if ~allHalfNotumIdenticalTF
    disp('ERROR in "RescaleAnimals": all animals listed in "avgAnimals" must correspond to same "halfNotum"!')
    return
end
% NB: this also ensures that "nMacroMAXspace" is the same for all animals.

% overwrites temp quantities
Animal = AnimalAOA;
PLOT = PLOTAOA;
temperature = temperatureAOA;



%% Loading archetype IF "Rescale" mode %%

if strcmp(rescaleMode,'Rescale')    
        
    % based on user input, select corresponding archetype
    if halfNotum == 'r'
        archName = 'HalfArchetype';
    elseif halfNotum == 'b'
        archName = 'FullArchetype';
    end
    
    files  = dir([archetypeFolder filesep archName '*.mat' ]);
    filename = files(1).name;
    archetype = load([archetypeFolder filesep files(1).name]);
    % NB: everything stored is in MICRONS
    
    uArchMacroX = archetype.macrochaete(1:end,1); % store X macrochaete coordinate (MICRONS)
    uArchMacroY = archetype.macrochaete(1:end,2); % store Y macrochaete coordinate (MICRONS)
    uArchMidPnt = archetype.midpoint;             % store Y midline coordinate (MICRONS)
    
    archA = [uArchMacroX uArchMacroY ; uArchMidPnt];
    archA(end,1) = NaN;
    stdArch   = [nanstd(uArchMacroX,0,2)       nanstd(uArchMacroY,0,2)];
    
elseif strcmp(rescaleMode,'Archetype') 
    archA = [];
end


%% Change of coordinate system: image (top left oginin) to archetype (new origin) %%

% Determine if mid landmarks coordinates (yMidline and xNeck) will be used
% in the determination of archetype and rescaling factors:
useMidPointTF = false;
%--------------------------------------------------------------------------
% % Use of the mid point into the origin calculation
% if all( ~isnan(dataNeck(3,:)) & ~isnan(dataMidLin(3,:)))
%   useMidPointTF = true;  
%   fprintf('warning: Mid landmarks will be use for origin calculation.\n');
% end
%--------------------------------------------------------------------------

[originX, originY] = ComputeOrigin(originType, dataMacroX, dataMacroY, dataMidLin, dataNeck, nMacroMin, halfNotum, useMidPointTF);
% NB: in PIXELS

% Coordinates refered to origin
oDataMacroX = dataMacroX - repmat(originX,[nMacroMAXspace 1]);  % PIXELS
oDataMacroY = dataMacroY - repmat(originY,[nMacroMAXspace 1]);  % PIXELS
oDataMidline = dataMidLin - repmat(originY,[3 1]);              % PIXELS
oDataNeck   = dataNeck   - repmat(originX,[3 1]);               % PIXELS


%% Determining rescaling factors %%

if ~exist(pathFolderRA,'dir')
    mkdir(pathFolderRA);
end

nMacroMAXspace = nMacroMAXspace + 1; % temporarily adds 1 to support (xNeck, yMidline) point

% Reshaping "scale1Dvec" into nMacroMAXspace x nAvgAnimals matrix:
scale1Dmat = repmat(scale1Dvec,1,nMacroMAXspace);
scale1Dmat = scale1Dmat';
% scale1Dmat = repmat(scale1D, nMacroMAXspace, nMovies);

Xall = [oDataMacroX ; oDataNeck(3,:)];   % add Neck coordinate (can be NaN)
Xall = Xall .* scale1Dmat;            % Convert PIXEL coordinates into MICRON

Yall = [oDataMacroY ; oDataMidline(3,:)]; % add Midline coordinate (can be NaN)
Yall = Yall .* scale1Dmat;            % Convert PIXEL coordinates into MICRON

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[moviesA, moviesB, scalingFactorX, scalingFactorY] = ComputeRescale(Xall, Yall, halfNotum, archA); % in MICRONS
% NB: if "archA" is empty, average of current animals will be used as archetype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update of "archA" that was empty if "Archetype" mode:
if strcmp(rescaleMode,'Archetype')
    archA = moviesA;
    archetype.range = zeros(4,1);
end

% Standard deviation around archetype position of Archetype and Movies to be rescale;
stdAnimals = [nanstd(Xall,0,2) nanstd(Yall,0,2)]; % 3.0
                         
% Standard deviation before and after rescaling; ** Use of Xall and Yall (1.3) **
stdNotRescaled = [ nanstd(Xall, 0, 2)  nanstd(Yall, 0, 2) ];
stdRescaled     = [ nanstd(Xall .* repmat(scalingFactorX,nMacroMAXspace,1),0,2) nanstd(Yall .* repmat(scalingFactorY,nMacroMAXspace,1),0,2) ];
                 
nMacroMAXspace = nMacroMAXspace - 1; % back to actual value (3.0)

                 
%% Save outputs

macroRawXY = [ dataMacroX ; dataMacroY ]; 

% rescaling output information
rescalingOutput = cell(8,nAvgAnimals+1);
rescalingOutput(1,1)       = {'Name'};
rescalingOutput(1,2:end)   = avgAnimals;
rescalingOutput(2:3,1)     = [{'Ox (pix)'}; {'Oy  (pix)'}];
rescalingOutput(2:3,2:end) = num2cell([originX ; originY]);
rescalingOutput(4:5,1)     = [{'xFactor'}; {'yFactor'}];
rescalingOutput(4:5,2:end) = num2cell([scalingFactorX ; scalingFactorY]);
rescalingOutput(6:7,1)     = [{'yML (pix)'}; {'xNK (pix)'}];
rescalingOutput(6:7,2:end) = num2cell([dataMidLin(3,:) ; dataNeck(3,:)]);
rescalingOutput(8,1)       = {'scale1D'};
rescalingOutput(8,2:end)   = num2cell(scale1Dmat(1,:));
% write output in .mat and .txt format
save([pathFolderRA filesep 'rescalingOutput' originTypeTag],'rescalingOutput');  
dlmcell([pathFolderRA filesep 'rescalingOutput' originTypeTag '.txt'],rescalingOutput,' ')

% std ouput information
stdOutput = cell(6,2);
stdOutput(1:2,1) = [{'stdAnimalsX (um)'} ; {'stdAnimalsY (um)'}];
stdOutput(1:2,2) = num2cell([nanmean(stdAnimals(:,1)) ; nanmean(stdAnimals(:,2))]);
stdOutput(3:4,1) = [{'stdNotRescaledX (um)'} ; {'stdNotRescaledY (um)'}];
stdOutput(3:4,2) = num2cell([nanmean(stdNotRescaled(:,1)) ; nanmean(stdNotRescaled(:,2))]);
stdOutput(5:6,1) = [{'stdRescaledX (um)'} ; {'stdRescaledY (um)'}];
stdOutput(5:6,2) = num2cell([nanmean(stdRescaled(:,1)) ; nanmean(stdRescaled(:,2))]);
% write output in .mat and .txt format
save([pathFolderRA filesep 'stdOutput' originTypeTag],'stdOutput');
dlmcell([pathFolderRA filesep 'stdOutput' originTypeTag '.txt'],stdOutput,' ')


%% Graphical Outputs

cmap = colormap(jet); close;
factor = size(cmap,1) / nAvgAnimals - 1;
size1 = 1:nMacroMAXspace;
size2 = [1:nMacroMAXspace 1:nMacroMAXspace];

% determination of max values for plot (1.2,1.3)
Xmin = min([min(Xall(:)) ; archetype.range(1)]);
Xmax = max([max(Xall(:)) ; archetype.range(2)]); % using "max" instead of "min" (4.0)
Ymax = max([max(abs(Yall(:))); archetype.range(4)]);
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
for n = 1:nAvgAnimals
    % color code
    color = cmap(round(factor * n),:);
    % scaling
    midlineVectorY = [oDataMidline(n) oDataMidline(n)] .* scale1Dmat(1,n); % (um)
    neckVectorX    = [oDataNeck(n) oDataNeck(n)]     .* scale1Dmat(1,n); % (um)
    macroVectorX   = oDataMacroX(:,n)                .* scale1Dmat(1,n); % (um)
    macroVectorY   = oDataMacroY(:,n)                .* scale1Dmat(1,n); % (um)
    % midline plot
    line([Xmin Xmax], midlineVectorY, 'LineStyle', ':', 'Color', color, 'LineWidth', 2)
    % macrochaete plot
    scatter(macroVectorX, macroVectorY, 20, color, 'LineWidth', 2);
    % neck plot
    scatter(neckVectorX, midlineVectorY, 20, color, 'LineWidth', 2);
end
% archetype plot
plot(archA(1:end-1,1), archA(1:end-1,2).', 'xb', 'LineWidth', 2, 'MarkerSize', 10)
plot(archA(end,1), archA(end,2), 'xb', 'LineWidth', 2, 'MarkerSize', 5)
line([Xmin Xmax], [archA(end,2) archA(end,2)], 'Color', [0 0 0], 'LineWidth', 2)
% origin (0,0)
plot(0, 0, '*k', 'LineWidth', 2, 'MarkerSize', 5)
% save figure
title('Raw Animals vs Archetype ({\mu}m)')
print('-dpng', '-r400', [pathFolderRA filesep 'raw_animals_vs_archetype' originTypeTag '.png']); % 1.1
if makeSVG
    plot2svg([pathFolderRA filesep 'raw_animals_vs_archetype' originTypeTag '.svg'], figure(1), 'png'); % 1.1
end
close

% Rescaled movies vs archetype
figure(2)
hold on
axis equal                  % 1.2
axis ij                     % 1.2
box on                      % 1.2
axis([Xmin Xmax Ymin Ymax]) % 1.2
for n = 1:nAvgAnimals
    % color code
    color = cmap(round(factor * n),:);
    % scaling and registering
    midlineVectorY = [oDataMidline(n) oDataMidline(n)] .* scalingFactorY(n) .* scale1Dmat(1,n); % (um)
    neckVectorX    = [oDataNeck(n)   oDataNeck(n)]   .* scalingFactorX(n) .* scale1Dmat(1,n); % (um)
    macroVectorX = scalingFactorX(n) .* oDataMacroX(:,n) .* scale1Dmat(1,n); % (um)
    macroVectorY = scalingFactorY(n) .* oDataMacroY(:,n) .* scale1Dmat(1,n); % (um)
    % mideline plot
    line([Xmin Xmax], midlineVectorY, 'LineStyle', ':', 'Color', color, 'LineWidth', 2)
    % macrochaete plot
    scatter(macroVectorX, macroVectorY, 20, color, 'LineWidth', 2)
    % neck plot
    scatter(neckVectorX, midlineVectorY, 20, color, 'LineWidth', 2);
end
% archetype plot
plot(archA(1:end-1,1), archA(1:end-1,2).', 'xb', 'LineWidth', 2, 'MarkerSize', 10)
plot(archA(end,1), archA(end,2), 'xb', 'LineWidth', 2, 'MarkerSize', 5)
line([Xmin Xmax], [archA(end,2) archA(end,2)], 'Color', [0 0 0], 'LineWidth', 2)
% origin (0,0)
plot(0, 0, '*k', 'LineWidth', 2, 'MarkerSize', 5)
% save figure
title('Rescaled Animals vs Archetype ({\mu}m)')
print('-dpng', '-r400', [pathFolderRA filesep 'rescaled_animals_vs_archetype' originTypeTag '.png']); % 1.1
if makeSVG
    plot2svg([pathFolderRA filesep 'rescale_animals_vs_archetype' originTypeTag '.svg'], figure(1), 'png'); % 1.1
end
close


%%  Saving txt file indicating date and version used in "Save_folder" (1.3)

today = datestr(now,29); % format 29 displays date yyyy-mm-dd style. Look up date for details
txtFilename = [today '_RescaleAnimals_' version '.txt'];
% Writing main parameters in txt file (3.1)
parameterCell = {   'Main Parameters:',[];
                    [],[];
                    'archetypeFolder = ', archetypeFolder;
                    'nMacroMin = ', nMacroMin;
                    'nAvgAnimals = ', nAvgAnimals;
                    'halfnotum = ', halfNotum;
                    'originType = ', originType;
                    'useMidPointTF = ', useMidPointTF};
                
dlmcell([pathFolderRA filesep txtFilename], parameterCell,' ');


%% Generate right and full archetype IF "Archetype" mode %%

if strcmp(rescaleMode,'Archetype')
    
    [Full, Half] = CreateSymmetricArchetype(moviesA(1:end-1,1),moviesA(1:end-1,2), moviesA(end,:), halfNotum, nMacroMin, useMidPointTF);
    % NB: everything is stored in MICRONS
    
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
    for a = 1:nAvgAnimals
        % we apply the shifted origin
        midlineVectorY = oDataMidline(a)  .* scale1Dmat(1,a) - Full.origin(2); % (um), mod 3.0
        neckVectorX    = oDataNeck(a)     .* scale1Dmat(1,a) - Full.origin(1); % (um), mod 3.0
        macroVectorX   = oDataMacroX(:,a) .* scale1Dmat(1,a) - Full.origin(1); % (um), mod 3.0
        macroVectorY   = oDataMacroY(:,a) .* scale1Dmat(1,a) - Full.origin(2); % (um), mod 3.0
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
    title('Full symmetric archetype ({\mu}m)')
    print('-dpng', '-r400', [pathFolderRA filesep 'full_symetric_archetype' originTypeTag '.png']); % 1.1
    if makeSVG
        plot2svg([pathFolderRA filesep 'full_symetric_archetype' originTypeTag '.svg'], figure(3), 'png'); % 1.1
    end
    close
    
    figure(3)
    hold on
    axis equal                  % 1.2
    axis ij                     % 1.2
    box on                      % 1.2
    axis([Xmin Xmax Ymin Ymax]) % 1.2
    for a = 1:nAvgAnimals
        % we apply the shifted origin
        midlineVectorY = oDataMidline(a)   .* scale1Dmat(1,a) - Half.origin(2); % (um)
        neckVectorX    = oDataNeck(a)     .* scale1Dmat(1,a) - Half.origin(1); % (um)
        macroVectorX   = oDataMacroX(:,a) .* scale1Dmat(1,a) - Half.origin(1); % (um)
        macroVectorY   = oDataMacroY(:,a) .* scale1Dmat(1,a) - Half.origin(2); % (um)
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
    title('Half symmetric archetype ({\mu}m)')
    print('-dpng', '-r400', [pathFolderRA filesep 'half_symetric_archetype' originTypeTag '.png']); % 1.1
    if makeSVG
        plot2svg([pathFolderRA filesep 'half_symetric_archetype' originTypeTag '.svg'], figure(3), 'png'); % 1.1
    end
    close
    
    % save archetype
    save([pathFolderRA filesep 'FullArchetype'],'-struct','Full');
    save([pathFolderRA filesep 'HalfArchetype'],'-struct','Half');
    % NB: everything is stored in MICRONS
    
end

%% History %%

% IMPROVEMENTS:
% - support cases where initial clicks starts with NaNs
% - support pupae oriented along vertical axis
% - add symmetry estimation during archetype calculation

% 25/04/2018: 4.0 (Boris) INCORPORATION OF "CreateArchetype"
% - accordingly defined additional parameter "rescaleMode" that either
% enable creation of archetype ("Archetype") or rescaling of animals 
% - now using "AIAFolderName" to find AIA_info files

% 25/04/2018: 3.0 (Boris)
% - now runs through "MAP_parameters"
% - now loads each animal "scale1D" to build "scale1Dmat" and therefore can
% support animal sets with different pixel scales.
% - now directly fetch txt files from landmarks clicks in each animal "SR"
% folders => no longer needs to gather them in "pathFolderArchetype" folder.
% - moved definition of "outputFolderCA" to "MAP_parameters"
% - removed "originType" in filenames, unless it is different than "2"
% (which we constantly use)

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

