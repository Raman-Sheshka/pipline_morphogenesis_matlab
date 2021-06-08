
Version = 3.6;
% Stephane Rigaud
% Anais Bailles
% Boris Guirao


%% Building of global common grid %%

% Merge all the GRID information loaded from the different annimals into
% the union of all the grids.

fprintf('\tBuilds global common grid ...')

% Load grids to define crop and translation
clear GRID;
nBoxes2OriginX = zeros(nAllAnimals, 1);
nBoxes2OriginY = zeros(nAllAnimals, 1);
lX = zeros(nAllAnimals, 1);
lY = zeros(nAllAnimals, 1);

for n = 1:nAllAnimals % for all the animal
    
    aGRID = AnimalGridList{n};
    
    if ~cloneTracking % 3.1
        
        % find the coordinates (integers) in box number of the barycenter (that normally defines the ULC of the [0,0] box in the grid)
        ULC = cell2mat3D(aGRID.ULCs);
        [~, Ox_a] = ind2sub(aGRID.Size, find( ULC(:,:,1) == aGRID.xywh(1)));
        [Oy_a, ~] = ind2sub(aGRID.Size, find( ULC(:,:,2) == aGRID.xywh(2)));
        nBoxes2OriginX(n) = unique(Ox_a); % box number (along x) of origin compartment (whose ULC = barycenter) starting from left side of grid
        nBoxes2OriginY(n) = unique(Oy_a); % box number (along y) of origin compartment (whose ULC = barycenter) starting from top side of grid
        % NB: nBoxes2OriginX/Y include the origin compartment (having [0,0] as "coordinates") in the count
        
        % find total number of compartments along x and y for this animal
        lX(n) = aGRID.Size(2);
        lY(n) = aGRID.Size(1);
        
    else % clone tracking case (3.2)
    
        % adding one row and one column all around data for display covenience (3.2)
        lX(n) = aGRID.Size(2);
        lY(n) = aGRID.Size(1);
        nBoxes2OriginX(n) = 1;
        nBoxes2OriginY(n) = 1;
        % NB: this averaging assumes that clone tracking was run with
        % "matchingWTclone" = 2, ONLY case where clone parts are gathered
        % together!
    end
end


% Find common global grid size
LtotX = max(nBoxes2OriginX) + max(lX - nBoxes2OriginX); % lx-Ox = number of boxes at the right of box Ox
LtotY = max(nBoxes2OriginY) + max(lY - nBoxes2OriginY); % ly-Oy = number of boxes below box Oy
gridSize = [LtotY LtotX];                             % manually defines "gridSize" to build it (2.7)

% Spot grid location in global common grid for each animal
startLocXwt = zeros(nAvgAnimals, 1);
startLocYwt = zeros(nAvgAnimals, 1);
endLocXwt = zeros(nAvgAnimals, 1);
endLocYwt = zeros(nAvgAnimals, 1);
for n = 1:nAvgAnimals
    [~,idx] = ismember(avgAnimals{n}, allAnimals);
    startLocXwt(n) = max(nBoxes2OriginX) - nBoxes2OriginX(idx) + 1; % +1 stems from the location being called from 1 to n (cf the piquets and barrieres pb)
    endLocXwt(n) = max(nBoxes2OriginX) - nBoxes2OriginX(idx) + lX(idx);
    startLocYwt(n) = max(nBoxes2OriginY) - nBoxes2OriginY(idx) + 1;
    endLocYwt(n) = max(nBoxes2OriginY) - nBoxes2OriginY(idx) + lY(idx);
end

% Determines image size that can contain the common grid, based on the box size:
% MINIMAL size of image that EXACTLY covers all grid compartments (2.7)
sizeImageXmin = ceil(DISPLAY.boxSize(1) * (1 + (LtotX-1) * (1 - gridOverlap)));     % added "ceil" (2.3), +1 added (2.7)
sizeImageYmin = ceil(DISPLAY.boxSize(2) * (1 + (LtotY-1) * (1 - gridOverlap)));     % added "ceil" and leeway of 5% (2.3), +1 added (2.7)
% NB: absolute minimum of sizeImage is the box size, even at overlap = 1 => 1+...
% NB: at olap = 0, image size min is exactly [LtotY*boxSize(2),LtotX*boxSize(1)]

%% Defining final image size to grid common grid (revamped 2.7) %%

% NB: +1 added so that, once we have set the origin that is NEVER (0,0), we STILL get gridSize = [LtotY LtotX]

% Adding at least one full box size along X and Y on BOTH sides of images (2.7)
leeway = 0.05; % 2.7
sizeImageXext = ceil(min(sizeImageXmin * (1+leeway), sizeImageXmin + 2 * DISPLAY.boxSize(1)));
sizeImageYext = ceil(min(sizeImageYmin * (1+leeway), sizeImageYmin + 2 * DISPLAY.boxSize(2)));
% NB: Adding leeway so as NOT to start grid at [0,0] but further from image border and still fit it in image

% setting final leeway:
leewayPixX = sizeImageXext - sizeImageXmin; % image extra pixels along X
leewayPixY = sizeImageYext - sizeImageYmin; % image extra pixels along Y

%% REG grid case (2.7) %%

sizeImageXreg = sizeImageXext;
sizeImageYreg = sizeImageYext;

% Setting new origin
% new origin with NO leeway (namely with grid starting at [0,0] of image)
newOxRaw = DISPLAY.boxSize(1) * (max(nBoxes2OriginX)-1) * (1 - gridOverlap);
newOyRaw = DISPLAY.boxSize(2) * (max(nBoxes2OriginY)-1) * (1 - gridOverlap);
% NB: subtracting "DISPLAY.boxSize" because box origin is normally at ULC of box => subtraction of one full box width/height
% ACTUAL new origin in image:
newOxReg = newOxRaw + leewayPixX/2;
newOyReg = newOyRaw + leewayPixY/2;
meanYmidReg = newOyReg - meanDeltaYmid; % 2.5

% Build global common REGULAR grid
% since using gridSize to keep previous size, need to determine grid new ULC in extended image (2.7)
gridULCx = newOxReg - DISPLAY.boxSize(1) * (max(nBoxes2OriginX) - 1) * (1 - gridOverlap);
gridULCy = newOyReg - DISPLAY.boxSize(2) * (max(nBoxes2OriginY) - 1) * (1 - gridOverlap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIDreg = MakeGrid([sizeImageYreg sizeImageXreg], DISPLAY.boxSize, [gridULCx gridULCy], gridSize , gridColor, gridLineWidth, gridOverlap);
% NB: at this stage, array "Coordinates" is NOT up to date as temporary
% grid origin was taken at top left corner of grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updating grid xywh by actual origin [newOxReg newOyReg]:
GRIDreg.xywh = [newOxReg newOyReg GRIDreg.xywh(3) GRIDreg.xywh(4)];                 
% Updating array "Coordinates" to set compartment containing orgin to [0,0] (3.4):
[gridCoordinates, originBoxIJ] = gridULCs2gridCoordinates(GRIDreg.ULCs, round([newOxReg newOyReg]));  
% NB: rounding [newOxReg newOyReg] to the nearest pixel to avoid assigning
% wrong compartement as origin.
GRIDreg.Coordinates = gridCoordinates;
GRIDreg.originBoxIJ = originBoxIJ;

% Macrochaetes:
macrocaetesReg = meanMacrocaetes; % in pixels, set for newly made up image
macrocaetesReg(1,:) = macrocaetesReg(1,:) + newOxReg;
macrocaetesReg(2,:) = macrocaetesReg(2,:) + newOyReg;

% Filling REG structure (2.7)
REG.xywh = GRIDreg.xywh;
REG.Centroids = GRIDreg.Centroids;
REG.Macrocaetes = macrocaetesReg;
REG.yMid = meanYmidReg;
REG.SizeImageX = sizeImageXreg;
REG.SizeImageY = sizeImageYreg;
% Saves macro & midline positions in MICRON with macrocaetes BARYCENTER AS ORIGIN (3.6)
REG.MacrocatesRaw = meanMacrocaetesRaw;
REG.yMidRaw = - meanDeltaYmidRaw;

%% Origin & Macrochaetes coordinates in MCU (Matrix Compartment Units) (3.5)%%

% NB: useful when directly plotting a matrix M using "image(M)" where 1
% compartment value is plotted with one pixel. Oddly, the first matrix
% compartment [1,1] is plotted between [0.5 0.5] and [1.5 1.5] (so that its
% center is at [1,1]), hence the [1/2 1/2] added at the origin.

newOxMCU = (originBoxIJ(2)-1) * (1 - gridOverlap) + 1/2;
newOyMCU = (originBoxIJ(1)-1) * (1 - gridOverlap) + 1/2;
macrocaetesMCU = meanMacrocaetes; % in pixels, set for newly made up image
macrocaetesMCU(1,:) = macrocaetesMCU(1,:)/DISPLAY.boxSize(1) + newOxMCU;
macrocaetesMCU(2,:) = macrocaetesMCU(2,:)/DISPLAY.boxSize(2) + newOyMCU;
% macrocaetesMCU = [[1;1] macrocaetesMCU]; % DEBUG

% Midline
yMidMCU = newOyMCU - meanDeltaYmid/DISPLAY.boxSize(2);

% Filling REG structure 
REG.MacrocaetesMCU = macrocaetesMCU;
REG.xywhMCU = [newOxMCU newOyMCU 1 1];
REG.yMidMCU = yMidMCU;


%% FULL grid case (2.7) %%

if ~strcmp(halfNotum,'b') % 2.8
    
    % Redefining newOy, sizeImageY, GRID and macrocaetes if making full animal out of halves (2.5)
    halfSizeImageYFull = sizeImageYreg - meanYmidReg;           % will only flip the part of image beyond midline (2.7)
    newOyFull = newOyReg + halfSizeImageYFull - meanYmidReg;    % 2.7
    newOxFull = newOxReg;                                       % keeps X coordinate (2.7)
    meanYmidFull = newOyFull - meanDeltaYmid;                   % meanYmid update: NOW RIGHT AT CENTER OF IMAGE ALONG Y (see NB below)
    % NB: meanYmidFull  = (newOyReg + halfSizeImageYFull - meanYmidReg) - meanDeltaYmid
    %                   = halfSizeImageYFull + newOyReg - newOyReg + meanDeltaYmid - meanDeltaYmid = halfSizeImageYFull
    sizeImageXFull = sizeImageXreg;                             % keeps image width (2.7)
    sizeImageYFull = ceil(2 * halfSizeImageYFull);              % 2.7
    
    % since using gridSize to keep previous size, need to determine grid new ULC in extended image (2.7)
    gridULCx = newOxFull - DISPLAY.boxSize(1) * (max(nBoxes2OriginX) - 1) * (1 - gridOverlap);
    gridULCy = newOyFull - DISPLAY.boxSize(1) * (max(nBoxes2OriginY) - 1) * (1 - gridOverlap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GRIDfull = MakeGrid([sizeImageYFull sizeImageXFull], DISPLAY.boxSize, [gridULCx gridULCy], gridSize , gridColor, gridLineWidth, gridOverlap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GRIDfull.xywh = [newOxFull newOyFull GRIDfull.xywh(3) GRIDfull.xywh(4)]; % replacing gridULCx/y by actual origin (adapted to full image)
    % Updating array "Coordinates" to set compartment containing orgin to [0,0] (3.4):
    [gridCoordinates, originBoxIJ] = gridULCs2gridCoordinates(GRIDreg.ULCs, round([newOxFull newOxFull]));
    % NB: rounding [newOxFull newOxFull] to the nearest pixel to avoid assigning
    % wrong compartement as origin.
    GRIDfull.Coordinates = gridCoordinates;
    GRIDfull.originBoxIJ = originBoxIJ;
    
    % Macrochaetes:
    macrocaetesFull = zeros(2,16);
    macrocaetesFull(1,1:8) = meanMacrocaetes(1,:) + newOxFull;
    macrocaetesFull(2,1:8) = meanMacrocaetes(2,:) + newOyFull;
    % filling coordinates for 8 flipped macrocaetes:
    macrocaetesFull(1,9:16) = macrocaetesFull(1,1:8);
    macrocaetesFull(2,9:16) = 2 * meanYmidFull - macrocaetesFull(2,1:8); % (yMacro + yMacroFlipped)/2 = meanYmid => yMacroFlipped = ...
    
    % Filling FULL structure (2.7)
    FULL.xywh = GRIDfull.xywh;
    FULL.Centroids = GRIDfull.Centroids;
    FULL.Macrocaetes = macrocaetesFull;
    FULL.yMid = meanYmidFull;
    FULL.SizeImageX = sizeImageXFull;
    FULL.SizeImageY = sizeImageYFull;
    
else
    FULL = REG; % when animals are "full" (halfNotum = 'b') (2.8)
end

%% Making a common GRID structure (2.7) %%

GRID = GRIDreg;
GRID.REG = REG;
GRID.FULL = FULL;
GRID = rmfield(GRID,{'xywh','Centroids'}); % those fields are in REG and FULL

fprintf(' Done\n')


% TOTO = cellfun(@(x)sum(abs(x - [newOxReg newOyReg])),GRID.ULCs,'UniformOutput',false);
% TATA = cell2mat(TOTO);

%% Initialization %%

% set up and test all variable for the averaging process
% NB: specific case if averaging AOS process for managing multichannel

fprintf('\tUniformizing data format ...')

temporary = AnimalAOTList{1}; % 2.3
list = fieldnames(temporary);
for m = 2:nAvgAnimals
    temporary = AnimalAOTList{m}; % 2.3
    list = intersect(list, fieldnames(temporary));
end
list = setdiff(list, rejectList);
nList = length(list); % 2.3

% dimension test
nbFrame = size(AnimalAOTList{1}.FrameArray, 1);% 2.3
testDim = 4;
if nbFrame == 1, testDim = 3; end

%%% Turning all 3D matrix (X,Y,Time) into 4D matrix (X,Y,1,Time)
for n = 1:nList
    for m = 1:nAvgAnimals
        eval( ['temporary = AnimalAOTList{m}.' list{n} ';' ] ); % 2.3
        
        if size(temporary,testDim) == 1
            new = reshape( temporary, size(temporary,1), size(temporary,2) ,1, size(temporary,3) );
            eval( ['AnimalAOTList{m}.' list{n} ' = new ;' ] ); % 2.3
        end
    end
end
fprintf(' Done\n')


%% Calculating averages over animals %%

% AOA backup initialisation
% We try to mimic the structure of single animal backup
BACKUP = GRID;
stdMAP = GRID;
fullBackup = GRID;

% find the number of time point
[Lt] = size( eval(['AnimalAOTList{1}.' list{1}]) , 4 ); % 2.3

%%% Stackin up quantities
% for each quantities, including AreaRatios and Rcons, we create a 5D-matrix (X,Y,Z,Time,Animal)
% wherw X and Y are defined by the new GRID, Z the dimension of the quantity values (1 scalar, 2 vector, or 4 tensor), Time the number of time point and Animal the number of animal
fprintf('\tStacking up data ...')
for m = 1:nList
    [~,~,Lz,~] = size( eval(['AnimalAOTList{1}.' list{m}]) ); % get the dimensions of the quantity (2.3)
    eval(['total' list{m} '= NaN(LtotY, LtotX, Lz, Lt, nAvgAnimals);']); % create a 5D matrix
    
    % Specific case for AreaRatios (2.3)
    if strcmp(list{m},'AreaRatios')
        totalAreaRatios = zeros(LtotY, LtotX, Lz, Lt, size(avgAnimals, 1)); % NB: RConds already included in AreaRatios for TA
    end
    
    for a = 1:nAvgAnimals
%         if ~cloneTracking
            eval(['total' list{m} '(startLocYwt(a):endLocYwt(a), startLocXwt(a):endLocXwt(a),:,:,a) = AnimalAOTList{a}.' list{m} ';']); % 2.3
%         else
%             eval(['total' list{m} '(startLocYwt(a)+1:endLocYwt(a)+1, startLocXwt(a)+1:endLocXwt(a)+1,:,:,a) = AnimalAOTList{a}.' list{m} ';']); % 3.2
%             eval(['total' list{m} '(2, 2,:,:,a) = AnimalAOTList{a}.' list{m} ';']); % ONLY fills compartment in the center, other compartments created just for display (3.1)
%         end
    end
end
fprintf(' Done\n')


%%% Calculating weights
fprintf('\tCalculating weights ...')
totalWeight = totalAreaRatios .^ 2;
% NB: RConds have been integrated to "total_AreaRatios" at "AIA_MultiOperation" stage
fprintf(' Done\n')


%%% Averaging values using weights when needed
fprintf('\tAveraging quantities ...\n')
for m = 1:nList
    fprintf(['\t\tProcessing ' list{m} ' ...'])
    
    if strcmp(list{m},'AreaRatios')
        meanTotalAreaRatios = mean(totalAreaRatios, 5);
        BACKUP.AreaRatios = meanTotalAreaRatios;
        stdMAP.AreaRatios = meanTotalAreaRatios;
        
        mergeBackup = cat(5, meanTotalAreaRatios, totalAreaRatios);
        eval(['fullBackup.' list{m} ' = mergeBackup;'])
    else
        
        currentTotalQuantity = eval(['total' list{m}]);
        % repeat weight matrix to fit the quantity structure (scalar, vector, or tensor)
        [Tz] = size(currentTotalQuantity, 3);
        R_TotalWeight = repmat(totalWeight, [1 1 Tz 1 1]);
        % calculate the weighted average of the quantity
        meanCurrentTotalQuantity = currentTotalQuantity .* R_TotalWeight;
        meanCurrentTotalQuantity = nansum( meanCurrentTotalQuantity, 5) ./ sum( R_TotalWeight, 5 );
        % calculate the significance map (binary) and std of the quantity
        [SignificanceMap, meanMap, stdMap] = MakeSignificanceMap( currentTotalQuantity, meanCurrentTotalQuantity, R_TotalWeight );
        if ~isempty(PLOT.SignOpacityMap)
            % calculate significance opacity value based on weights
            SignOpacityMap_iso = max(SignificanceMap(:,:,1,:) , abs(SignificanceMap(:,:,1,:) - 1) .* PLOT.SignOpacityMap(1) );
            SignOpacityMap_dev = max(SignificanceMap(:,:,2,:) , abs(SignificanceMap(:,:,2,:) - 1) .* PLOT.SignOpacityMap(2) );
            SignificanceMap = cat(3, SignificanceMap, SignOpacityMap_iso, SignOpacityMap_dev);
        end
        % store mean in new AOA backup
        eval(['BACKUP.' list{m} ' = meanCurrentTotalQuantity;']);
        eval(['BACKUP.' list{m} '_Smap = SignificanceMap;']);
        eval(['stdMAP.' list{m} '_std = meanMap;']);  % (v2) in progress add-on for biological variability estimation
        eval(['stdMAP.' list{m} ' = meanCurrentTotalQuantity;']);   % (v2) in progress add-on for biological variability estimation
        
        mergeBackup = cat(5, meanCurrentTotalQuantity, currentTotalQuantity);
        eval(['fullBackup.' list{m} ' = mergeBackup;'])
        eval(['fullBackup.' list{m} '_Smap = SignificanceMap;']);
    end
    fprintf(' Done\n')
end % end for each quantities


%%% Add the rest to the backup (mod 2.3)
BACKUP.TimeArray  = AnimalAOTList{1}.TimeArray;
BACKUP.FrameArray = AnimalAOTList{1}.FrameArray;
stdMAP.TimeArray  = AnimalAOTList{1}.TimeArray;  % (v2) in progress add-on for biological variability estimation
stdMAP.FrameArray = AnimalAOTList{1}.FrameArray; % (v2) in progress add-on for biological variability estimation

fullBackup.TimeArray = AnimalAOTList{1}.TimeArray;
fullBackup.FrameArray = AnimalAOTList{1}.FrameArray;

%%% Save AOA mean backup
if ~exist(OutputPathName,'dir')
    mkdir(OutputPathName);
end

BACKUP = fullBackup;

%%% Save main backup
thisFilename = [OutputPathName filesep 'AOA_backup' '.mat']; % added fullTag (2.6)
if ~exist(thisFilename,'file')
    fprintf('\tSaving new backup ...')
    save(thisFilename, '-struct', 'BACKUP');
    % Backup being tested
    %     save([OutputPathName filesep 'AOA_stdMAP' '.mat'], '-struct', 'stdMAP'); % added fullTag (2.6)
    %     save([OutputPathName filesep 'AOA_fullBackup' '.mat'], '-struct', 'fullBackup'); % added fullTag (2.6)
else
    fprintf('\tAppending backup ...')
    Previous = load(thisFilename);
    BACKUP = catstruct(Previous, BACKUP);
    save(thisFilename, '-struct', 'BACKUP', '-append');
    % Backup being tested
    %     save([OutputPathName filesep 'AOA_stdMAP' '.mat'], '-struct', 'stdMAP', '-append'); % added fullTag (2.6)
    %     save([OutputPathName filesep 'AOA_fullBackup' '.mat'], '-struct', 'fullBackup', '-append'); % added fullTag (2.6)
end
fprintf(' Done\n')

%% Plot %%

if ~isempty(Qname) && PLOT.plot && ~skipAOAplots
    fprintf('Plotting averaged quantities ...\n')
    
    BACKUP = load([OutputPathName filesep 'AOA_backup']); % added fullTag (2.6)
    % BACKUP = load([OutputPathName filesep 'AOA_backup_MOD']);
    % NB: "MOD" backups were made to symmetrize tensors ONLY in the row right next to the midline,
    % to mimick the calculation done in BIGwt2 were a big box (similar to Jesus' box) including
    % both sides of the midline can be drawn.
    
    %%% Loading plot parameters (mod 2.9)
    % Extracting FULL or REG data for plot (2.7)
    TYPE = BACKUP.REG;
    if DISPLAY.makeItFull
        TYPE = BACKUP.FULL;
    end
    BACKUP.xywh = TYPE.xywh;
    BACKUP.Centroids = TYPE.Centroids;
    BACKUP.SizeImageX = TYPE.SizeImageX;
    BACKUP.SizeImageY = TYPE.SizeImageY;
    DISPLAY.yMid = TYPE.yMid;
    DISPLAY.macrocaetes = TYPE.Macrocaetes;
    
    Qimage = ones(BACKUP.SizeImageY, BACKUP.SizeImageX);
    
    %%% For each time point
    for t = 1:Lt
        
        % Correcting time display: ONLY relevant when TA quantities were added LAST in AOT backup(2.5)
        fixDtH = 0;
        timeShift = 0;
        % fixDtH = delta_t/60; % uncomment if 12h05 is displayed instead of 12h00 for instance
        
        % Set the time and n variable (mod 2.2)
        tmpTimeStart = TimeStr2Dec(BACKUP.TimeArray{t,1}) + timeShift - fixDtH;      % turns time string into decimal and adds timeShift, added fixDtH (2.5)
        tmpTimeEnd = TimeStr2Dec(BACKUP.TimeArray{t,2}) + timeShift;
        tmpTimeStart = TimeDec2Str(tmpTimeStart);                                            % switches back to string
        tmpTimeEnd = TimeDec2Str(tmpTimeEnd);
        DISPLAY.time = [tmpTimeStart ' - ' tmpTimeEnd];
        DISPLAY.n = round( ( BACKUP.FrameArray(t,1) + BACKUP.FrameArray(t,2) ) / 2 );
        
        if Lt > 1 % (1.2)
            DISPLAY.step = t;
        end
        
        
        %%% For each Display to be plot
        for j = 1:nQname
            
            if isfield(BACKUP,Qname{j}) % 2.3
                
                % Get the Qcolor, Qunits, Qsr, etc. corresponding to the Display
                [Pname,idx] = GetPname(Qname{j}); % now gets Pname (2.3)
                Qunits    = eval(['allUnits' Pname '{idx}']);
                Qcolor    = eval(['allColors' Pname '{idx}']);
                Qsr       = eval(['scaleRatio' Pname ]);
                Qscalebar = eval(['scaleBarLength' Pname ]);
                QKillTr   = eval(['killMeanTrace' Pname ]);
                Qsr = GetPlotScaleParameter(Qsr, idx);
                Qscalebar = GetPlotScaleParameter(Qscalebar, idx);
                QKillTr = GetPlotScaleParameter(QKillTr, idx);
                
                % Get the Qcolor, Qunits, Qsr, etc. corresponding to the Display
                [Pname,idx] = GetPname(Qname{j}); % now gets Pname (2.3)
                Qunits    = eval(['allUnits' Pname '{idx}']);
                Qcolor    = eval(['allColors' Pname '{idx}']);
                
                
                
                % tag indicating mean Tr = 0 for naming files
                KillTrTag = '';
                if QKillTr
                    KillTrTag = '_Tr=0';
                end
                
                % Load significance map
                if PLOT.significance && ~singleAnimal
                    DISPLAY.SignificanceMap = eval(['BACKUP.' Qname{j} '_Smap(:,:,:,t)']);
                end
                
                % Plot image
                PlotField(Qname{j}, QKillTr, Qcolor , Qunits , Qsr, Qscalebar, BACKUP, Qimage, DISPLAY);
                
                % Plot print (removed "Frames_" prefix and added fullTag in filenames)(2.6)
                thisFilename = ['AOA_' Qname{j} '_' tmpTimeStart 'to' tmpTimeEnd KillTrTag fullTag '_sr=' num2str(Qsr(1)) imageExtension]; % 1.7, added fullTag (2.5), use timeStart/End (2.7)
                plotFolder = [OutputPathName filesep 'AOA_' plotType fullTag];
                if ~exist(plotFolder, 'dir')
                    mkdir(plotFolder);
                end
                if strcmp(imageExtension,'.svg')
                    plot2svg([plotFolder filesep thisFilename],figure(1),'png');
                else
                    print(printFormat, printResolution, [plotFolder filesep thisFilename]);
                end
                close
                
            else % 2.3
                disp(['AOA WARNING: quantity "' Qname{j} '" was not found in "BACKUP"!']);
                disp('Please make sure that it is the name assigned to this quantity AND that it was averaged over time and stored in each animal AOT backups.')
            end
        end % end for each quantities
    end % end for each time point
elseif ~isempty(Qname) && PLOT.plot && skipAOAplots
    disp('WARNING: skipping AOA plots ("skipAOAplots" set to 1)!')
end % end if plot


%% History %%

% 29/09/2020: 3.6 (Boris)
% - Saves macro & midline positions in MICRON with macrocaetes BARYCENTER AS ORIGIN!
% Saved in REG substructure (REG.MacrocatesRaw, REG.yMidRaw).

% 09/07/2020: 3.5 (Boris)
% - added the calculation (and storage in REG structure) of MacrocaetesMCU,
% xywhMCU and yMidMCU, coordinates of the corresponding quantities in MCU
% (matrix compartement units) for Lorette and Eric's plots.

% 02/06/2020: 3.4 (Boris)
% - call to " SignificanceMap" rather than "SignificantMap"

% 20/04/2020: 3.3 (Boris)
% - call to " SignificanceMap" rather than "SignificavityMap"

% 06/07/2018: 3.2 (Boris)
% - BETA: alternative averaging of clone that assumes that clone tracking
% has been done with "matchingWTclone" = 2, only case that gathers all
% clone parts together in the average.

% 02/07/2018: 3.1 (Boris)
% - supports clone averaging (without WT counterpart) and display.

% 13/06/2018: 3.0 (Stephane)
% - update for new SAP and MAP script
% - add a fifth dimension to compute average animal (A1) and to keep Qs for each animal (A2 to An)
% - remove single animal processing (single animal only done in SAP now)

% 08/02/2017: 2.9 (Boris)
% - Adjustments for single movie processing
% - added "skipAOAplots"

% 03/02/2017: 2.8 (Boris)
% - do NOT try to build full grid when processing animals tha are halfNotum = 'b'

% 23-26/01/2017: 2.7 (Boris)
% - now will save both regular and full versions of xywh, yMid, sizeImageX/Y, macrochaetes in stuctures REG and FULL within GRID that get
% extracted only right before plot => back to ONLY saving one backup for AOA
% - fixed the issue of mismatch between grid size and actual matrix size by defining gridSize right from LtotX/Y
% - fixed all formulas wrongly setting origins and sizes due to wrong handling of the overlap between boxes.
% - impose a minimum distance of either one full box everywhere around data or 5%/2 of image size (whichever is smallest)
% - now takes into account "timeShift" when naming files being saved

% 20/01/2017: 2.6 (Boris)
% - added a 5% leeway in sizeImageX to start drawing tensors further from left border
% - added fullTag to folders and backup names to make follwing programs compatible
% - took "fullTag" definition to AIA_MultiOperation
% - moved most filling of DISPLAY to AOA_MultiOperation

% 18/01/2017: 2.5 (Boris)
% - made plot of FULL mean animal out of a HALF mean animal possible
% - removed leeway along Y to define "sizeImageY" (issue before was due to something else)
% - removed plot of miline that is now done in PlotField

% 18/01/2017: 2.4 (Boris)
% - same as 2.3 but removed most of the many old commented parts
% - decreased leeway along Y to define "sizeImageY" to 5% (from 10%)

% 12/01/2017: 2.3 (Boris)
% - ONLY using AOT backups now!! (=> removed "PnameAOA")
% - added leeway of 10% along Y to define "sizeImageY"

% 05/01/2016: 2.2 (Boris)
% - now uses "timeShift" to correct displayed time range of averaging
% - "All_Animals" became "allAnimals"
% - "Animals" became "avgAnimals"
% - "genericName" became "Animal"
% - "Pname" became "PnameAOA"

% 26/05/2016: 2.1
% - RConds are now mean using nanmean to avoid them impacting weight when they are not calculated
%   This follow the AOT modification to allow time windows average regardless of the movie length
% - Verification after the weight calculation to switch NaN to 0, because weight shoudl alway be between {0;1}

% 25/03/2015: 2.0
% - reshape of the mean process for better clarity and efficiency. It now in two times:
%   1 - create for each quantities a 5D matrix X,Y,Z,Time,Animal
%   2 - apply mean on the 5th dimension of each matrix
% - include of RConds variable in TA process

% 21/01/2015: 1.2
% - fixed bug when nb time point data == 1

% 13/11/2014 v1.1 :
% - add the new plot_field procedure that allow to plot iso and aniso value at different scale
% - start of integration of VMM process. TC,SO,DRP are not tested
% - modification of the code to better fit AIA_parameters from Boris

% 22/10/2014 v1 :
% - simplified rip-off of Animal_Averager of Anais adapted to Boris framework