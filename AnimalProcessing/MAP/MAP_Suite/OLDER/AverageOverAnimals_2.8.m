% AverageOverAnimals (AOA)
%
% Script containing process to average data from multiple AOT runs
%
% Version 2.8
% Stephane Rigaud
% Anais Bailles
% Boris Guirao


%% Building of global common grid %%

% Merge all the GRID information loaded from the different annimals in mean_GRID_specs
% Needed for determining the grid on which to mean and the plot are made.

disp('Builds global common grid ...')

% Load grids to define crop and translation
clear GRID;
nBoxes2OriginX = zeros(nAllAnimals,1);
nBoxes2OriginY = zeros(nAllAnimals,1); 
lX = zeros(nAllAnimals,1); 
lY = zeros(nAllAnimals,1); 

for n = 1:nAllAnimals % for all the animal
    grid = mean_GRID_specs{n};
    
    % find the coordinates (integers) in box number of the barycenter (that normally defines the ULC of the [0,0] box in the grid)
    ULC = cell2mat3D(grid.ULCs);
    [useless1 Ox_a] = ind2sub(grid.size, find( ULC(:,:,1) == grid.xywh(1)));
    [Oy_a useless2] = ind2sub(grid.size, find( ULC(:,:,2) == grid.xywh(2)));
    nBoxes2OriginX(n) = unique(Ox_a); % box number (along x) of origin compartment (whose ULC = barycenter) starting from left side of grid
    nBoxes2OriginY(n) = unique(Oy_a); % box number (along y) of origin compartment (whose ULC = barycenter) starting from top side of grid
    % NB: nBoxes2OriginX/Y include the origin compartment (having [0,0] as "coordinates") in the count
    
    % find total number of compartments along x and y for this animal
    lX(n) = grid.size(2);
    lY(n) = grid.size(1);
end

% Find common global grid size
LtotX = max(nBoxes2OriginX) + max(lX-nBoxes2OriginX); % lx-Ox = number of boxes at the right of box Ox
LtotY = max(nBoxes2OriginY) + max(lY-nBoxes2OriginY); % ly-Oy = number of boxes below box Oy
gridSize = [LtotY LtotX];                             % manually defines "gridSize" to build it (2.7)

% Spot grid location in global common grid for each animal
startLocXwt = zeros(nAvgAnimals,1);
startLocYwt = zeros(nAvgAnimals,1);
endLocXwt = zeros(nAvgAnimals,1);
endLocYwt = zeros(nAvgAnimals,1);
for n = 1:nAvgAnimals
    [~,idx] = ismember(avgAnimals{n}, allAnimals);
    startLocXwt(n) = max(nBoxes2OriginX) - nBoxes2OriginX(idx) + 1; % +1 stems from the location being called from 1 to n (cf the piquets and barrieres pb)
    endLocXwt(n) = max(nBoxes2OriginX) - nBoxes2OriginX(idx) + lX(idx);
    startLocYwt(n) = max(nBoxes2OriginY) - nBoxes2OriginY(idx) + 1;
    endLocYwt(n) = max(nBoxes2OriginY) - nBoxes2OriginY(idx) + lY(idx);
end

% Determines image size that can contain the common grid, based on the box size:
% MINIMAL size of image that EXACTLY covers all grid compartments (2.7)
sizeImageXmin = ceil(PLOT.boxSize(1)*(1 + (LtotX-1)*(1 - gridOverlap)));     % added "ceil" (2.3), +1 added (2.7)
sizeImageYmin = ceil(PLOT.boxSize(2)*(1 + (LtotY-1)*(1 - gridOverlap)));     % added "ceil" and leeway of 5% (2.3), +1 added (2.7)
% NB: absolute minimum of sizeImage is the box size, even at overlap = 1 => 1+...
% NB: at olap = 0, image size min is exactly [LtotY*boxSize(2),LtotX*boxSize(1)]


%% Defining final image size to grid common grid (revamped 2.7) %%

% NB: +1 added so that, once we have set the origin that is NEVER (0,0), we STILL get gridSize = [LtotY LtotX]

% Adding at least one full box size along X and Y on BOTH sides of images (2.7)
leeway = 0.05; % 2.7
sizeImageXext = ceil(min(sizeImageXmin*(1+leeway), sizeImageXmin + 2*PLOT.boxSize(1)));
sizeImageYext = ceil(min(sizeImageYmin*(1+leeway), sizeImageYmin + 2*PLOT.boxSize(2)));
% NB: Adding leeway so as NOT to start grid at [0,0] but further from image border and still fit it in image

% sizeImageXext = ceil(sizeImageXmin)+1; % DEBUG
% sizeImageYext = ceil(sizeImageYmin)+1; % DEBUG

% setting final leeway:
leewayPixX = sizeImageXext - sizeImageXmin; % image extra pixels along X
leewayPixY = sizeImageYext - sizeImageYmin; % image extra pixels along Y


%% REG grid case (2.7) %%

sizeImageXreg = sizeImageXext;
sizeImageYreg = sizeImageYext;

% Setting new origin
% new origin with NO leeway (namely with grid starting at [0,0] of image)
newOxRaw = PLOT.boxSize(1) * (1 + (max(nBoxes2OriginX)-1) * (1 - gridOverlap)) - PLOT.boxSize(1);
newOyRaw = PLOT.boxSize(2) * (1 + (max(nBoxes2OriginY)-1) * (1 - gridOverlap)) - PLOT.boxSize(2);
% NB: subtracting "PLOT.boxSize" because box origin is normally at ULC of box => subtraction of one full box width/height
% ACTUAL new origin in image:
newOxReg = newOxRaw + leewayPixX/2;
newOyReg = newOyRaw + leewayPixY/2;
meanYmidReg = newOyReg - meanDeltaYmid; % 2.5

% Build global common REGULAR grid
% since using gridSize to keep previous size, need to determine grid new ULC in extended image (2.7)
gridULCx = newOxReg - PLOT.boxSize(1) * (max(nBoxes2OriginX)-1) * (1 - gridOverlap);
gridULCy = newOyReg - PLOT.boxSize(2) * (max(nBoxes2OriginY)-1) * (1 - gridOverlap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIDreg = GridMaker([sizeImageYreg sizeImageXreg], PLOT.boxSize, [gridULCx gridULCy], gridSize , gridColor, gridLineWidth, gridOverlap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIDreg.xywh = [newOxReg newOyReg GRIDreg.xywh(3) GRIDreg.xywh(4)]; % replacing gridULCx/y by actual origin

% Macrochaetes:
macrocaetesReg = DISPLAY.macrocaetes;
macrocaetesReg(1,:) = macrocaetesReg(1,:) + newOxReg;
macrocaetesReg(2,:) = macrocaetesReg(2,:) + newOyReg;

% Filling REG structure (2.7)
REG.xywh = GRIDreg.xywh;
REG.centroids = GRIDreg.centroids;
REG.macrocaetes = macrocaetesReg;
REG.yMid = meanYmidReg;
REG.sizeImageX = sizeImageXreg;
REG.sizeImageY = sizeImageYreg;


%% FULL grid case (2.7) %%

if ~strcmp(RAW.halfNotum,'b') % 2.8

    % Redefining newOy, sizeImageY, GRID and macrocaetes if making full animal out of halves (2.5)
    halfSizeImageYFull = sizeImageYreg - meanYmidReg;           % will only flip the part of image beyond midline (2.7)
    newOyFull = newOyReg + halfSizeImageYFull - meanYmidReg;    % 2.7
    newOxFull = newOxReg;                                       % keeps X coordinate (2.7)
    meanYmidFull = newOyFull - meanDeltaYmid;                   % meanYmid update: NOW RIGHT AT CENTER OF IMAGE ALONG Y (see NB below)
    % NB: meanYmidFull  = (newOyReg + halfSizeImageYFull - meanYmidReg) - meanDeltaYmid
    %                   = halfSizeImageYFull + newOyReg - newOyReg + meanDeltaYmid - meanDeltaYmid = halfSizeImageYFull
    sizeImageXFull = sizeImageXreg;                             % keeps image width (2.7)
    sizeImageYFull = ceil(2*halfSizeImageYFull);                % 2.7
    
    % since using gridSize to keep previous size, need to determine grid new ULC in extended image (2.7)
    gridULCx = newOxFull - PLOT.boxSize(1) * (max(nBoxes2OriginX)-1) * (1 - gridOverlap);
    gridULCy = newOyFull - PLOT.boxSize(1) * (max(nBoxes2OriginY)-1) * (1 - gridOverlap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GRIDfull = GridMaker([sizeImageYFull sizeImageXFull], PLOT.boxSize, [gridULCx gridULCy], gridSize , gridColor, gridLineWidth, gridOverlap);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GRIDfull.xywh = [newOxFull newOyFull GRIDfull.xywh(3) GRIDfull.xywh(4)]; % replacing gridULCx/y by actual origin (adapted to full image)
    
    % Macrochaetes:
    macrocaetesFull = DISPLAY.macrocaetes;
    macrocaetesFull(1,:) = macrocaetesFull(1,:) + newOxFull;
    macrocaetesFull(2,:) = macrocaetesFull(2,:) + newOyFull;
    % filling coordinates for 8 flipped macrocaetes:
    macrocaetesFull(1,9:16) = macrocaetesFull(1,1:8);
    macrocaetesFull(2,9:16) = 2*meanYmidFull - macrocaetesFull(2,1:8); % (yMacro + yMacroFlipped)/2 = meanYmid => yMacroFlipped = ...
    
    % Filling FULL structure (2.7)
    FULL.xywh = GRIDfull.xywh;
    FULL.centroids = GRIDfull.centroids;
    FULL.macrocaetes = macrocaetesFull;
    FULL.yMid = meanYmidFull;
    FULL.sizeImageX = sizeImageXFull;
    FULL.sizeImageY = sizeImageYFull;

else
    FULL = REG; % when animals are "full" (halfNotum = 'b') (2.8)
end

%% Making a common GRID structure (2.7) %%

GRID = GRIDreg;
GRID.REG = REG;
GRID.FULL = FULL;
GRID = rmfield(GRID,{'xywh','centroids'}); % those fields are in REG and FULL


%% Initialization %%

% set up and test all variable for the averaging process
% NB: specific case if averaging AOS process for managin multichannel

fprintf('Uniformizing data format ...')

eval( ['temporary = GRID_AOT_' avgAnimals{1} ';' ] ); % 2.3
list = fieldnames(temporary);
for m = 2:nAvgAnimals
    eval( ['temporary = GRID_AOT_' avgAnimals{m} ';' ] ); % 2.3
    list = intersect(list,fieldnames(temporary));
end
list = setdiff(list,rejectList);
nList = length(list); % 2.3

% dimension test
eval( ['nbFrame = size(GRID_AOT_' avgAnimals{1} '.FrameArray,1 );' ] );% 2.3
testDim = 4;
if nbFrame == 1, testDim = 3; end

%%% Turning all 3D matrix (X,Y,Time) into 4D matrix (X,Y,1,Time)
for n = 1:nList
    for m = 1:nAvgAnimals
        eval( ['temporary = GRID_AOT_' avgAnimals{m} '.' list{n} ';' ] ); % 2.3
        
        if size(temporary,testDim) == 1
            new = reshape( temporary, size(temporary,1), size(temporary,2) ,1, size(temporary,3) );
            eval( ['GRID_AOT_' avgAnimals{m} '.' list{n} ' = new ;' ] ); % 2.3
        end
    end
end
fprintf('Done\n')

disp('FYI: Loaded backup to be averaged ...')
for m = 1:nAvgAnimals
    disp(backupPathList{m}); % display the backup being loaded
end


%% Calculating averages over animals %%

% AOA backup initialisation
% We try to mimic the structure of single animal backup
AOA_backup = GRID;
AOA_stdMAP = GRID;

% find the number of time point
[Lt] = size( eval(['GRID_AOT_' avgAnimals{1} '.' list{1}]) , 4 ); % 2.3

%%% Stackin up quantities
%----------------------------------------------------------------------------------------
% for each quantities, including AreaRatios and Rcons, we create a 5D-matrix (X,Y,Z,Time,Animal)
% wherw X and Y are defined by the new GRID, Z the dimension of the quantity values (1 scalar, 2 vector, or 4 tensor), Time the number of time point and Animal the number of animal
fprintf('Stacking up data ...')
for m = 1:nList
    [~,~,Lz,~] = size( eval(['GRID_AOT_' avgAnimals{1} '.' list{m}]) ); % get the dimensions of the quantity (2.3)
    eval(['total_' list{m} '= NaN(LtotY, LtotX, Lz, Lt, nAvgAnimals);']); % create a 5D matrix
    
    % Specific case for AreaRatios (2.3)
    if strcmp(list{m},'AreaRatios')
        total_AreaRatios = zeros(LtotY, LtotX, Lz, Lt, size(avgAnimals,1)); % NB: RConds already included in AreaRatios for TA
    end
    
    for a = 1:nAvgAnimals
        eval(['total_' list{m} '(startLocYwt(a):endLocYwt(a), startLocXwt(a):endLocXwt(a),:,:,a) = GRID_AOT_' avgAnimals{a} '.' list{m} ';']); % 2.3
    end
end
fprintf('Done\n')
%----------------------------------------------------------------------------------------

%%% Calculating weights
%----------------------------------------------------------------------------------------
fprintf('Calculating weights ...')
total_weight = total_AreaRatios .^ 2;
% NB: RConds have been integrated to "total_AreaRatios" at "AIA_MultiOperation" stage
fprintf('Done\n')
%----------------------------------------------------------------------------------------

%%% Averaging values using weights when needed
%----------------------------------------------------------------------------------------
disp('Averaging quantities...')
for m = 1:nList
    
    if strcmp(list{m},'AreaRatios')
        disp(['    Processing ' list{m} '...']);
       AOA_backup.AreaRatios = mean(total_AreaRatios, 5);
       AOA_stdMAP.AreaRatios = mean(total_AreaRatios, 5);
         
    else
        disp(['    Processing ' list{m} '...']);
        % repeat weight matrix to fit the quantity structure (scalar, vector, or tensor)
        [Tz] = eval(['size( total_' list{m} ',3);']);
        R_total_weight = repmat(total_weight, [1 1 Tz 1 1]);
        % calculate the weighted average of the quantity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp_mean = eval(['total_' list{m} ' .* R_total_weight;']);
        temp_mean = nansum( temp_mean, 5) ./ sum( R_total_weight, 5 );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate the significance map (binary) and std of the quantity
        [SignificantMap,meanMap,stdMap] = SignificativityMap( eval(['total_' list{m}]), temp_mean, R_total_weight );
        if ~isempty(PLOT.SignOpacityMap)
            % calculate significance opacity value based on weights
            SignOpacityMap_iso = max(SignificantMap(:,:,1,:) , abs(SignificantMap(:,:,1,:) - 1) .* PLOT.SignOpacityMap(1) );
            SignOpacityMap_dev = max(SignificantMap(:,:,2,:) , abs(SignificantMap(:,:,2,:) - 1) .* PLOT.SignOpacityMap(2) );
            SignificantMap = cat(3, SignificantMap, SignOpacityMap_iso, SignOpacityMap_dev);
        end
        % store mean in new AOA backup
        eval(['AOA_backup.' list{m} ' = temp_mean;']);
        eval(['AOA_backup.' list{m} '_Smap = SignificantMap;']);
        eval(['AOA_stdMAP.' list{m} '_std = meanMap;']);  % (v2) in progress add-on for biological variability estimation
        eval(['AOA_stdMAP.' list{m} ' = temp_mean;']);   % (v2) in progress add-on for biological variability estimation
        
    end
end % end for each quantities
fprintf('Done\n')
%----------------------------------------------------------------------------------------

%%% Add the rest to the backup (mod 2.3)
eval(['AOA_backup.TimeArray  = GRID_AOT_' avgAnimals{1} '.TimeArray;']);
eval(['AOA_backup.FrameArray = GRID_AOT_' avgAnimals{1} '.FrameArray;']);
eval(['AOA_stdMAP.TimeArray  = GRID_AOT_' avgAnimals{1} '.TimeArray;']);  % (v2) in progress add-on for biological variability estimation
eval(['AOA_stdMAP.FrameArray = GRID_AOT_' avgAnimals{1} '.FrameArray;']); % (v2) in progress add-on for biological variability estimation

%%% Save AOA mean backup
if ~exist(OutputPathName,'dir')
    mkdir(OutputPathName); 
end
save([OutputPathName filesep 'AOA_stdMAP' '.mat'], 'AOA_stdMAP'); % added fullTag (2.6)

thisFilename = [OutputPathName filesep 'AOA_backup' '.mat']; % added fullTag (2.6)
if ~exist(thisFilename,'file')
    disp('Saving new backup ...')
    save(thisFilename, '-struct', 'AOA_backup');
else
    disp('Appending backup ...')
    Previous = load(thisFilename);
    AOA_backup = catstruct(Previous,AOA_backup);
    save(thisFilename, '-struct', 'AOA_backup', '-append');
end


%% Plot %%

if ~isempty(Qname) && PLOT.plot
    
    disp('Starting plot ...')
    
    AOA_backup = load([OutputPathName filesep 'AOA_backup']); % added fullTag (2.6)
%     AOA_backup = load([OutputPathName filesep 'AOA_backup_MOD']);
    % NB: "MOD" backups were made to symmetrize tensors ONLY in the row right next to the midline, 
    % to mimick the calculation done in BIGwt2 were a big box (similar to Jesus' box) including
    % both sides of the midline can be drawn.
    
    % Extracting FULL or REG data for plot (2.7)
    TYPE = REG;
    if PLOT.makeItFull
        TYPE = FULL;
    end
    AOA_backup.xywh = TYPE.xywh;
    AOA_backup.centroids = TYPE.centroids;
    AOA_backup.sizeImageX = TYPE.sizeImageX;
    AOA_backup.sizeImageY = TYPE.sizeImageY;
    
    DISPLAY.yMid = TYPE.yMid;
    DISPLAY.macrocaetes = TYPE.macrocaetes;
    
    Qimage = ones(AOA_backup.sizeImageY, AOA_backup.sizeImageX);

    %%% For each time point
    for i = 1:Lt
       
        % Correcting time display: ONLY relevant when TA quantities were added LAST in AOT backup(2.5)
        fixDtH = 0; 
%         fixDtH = delta_t/60; % uncomment if 12h05 is displayed instead of 12h00 for instance
        
        % Set the time and n variable (mod 2.2)
        timeStart = Time_str2dec(AOA_backup.TimeArray{i,1}) + timeShift - fixDtH;      % turns time string into decimal and adds timeShift, added fixDtH (2.5)
        timeEnd = Time_str2dec(AOA_backup.TimeArray{i,2}) + timeShift;
        timeStart = Time_dec2str(timeStart);                                            % switches back to string
        timeEnd = Time_dec2str(timeEnd);
        DISPLAY.time = [timeStart ' - ' timeEnd];
        DISPLAY.n = round( ( AOA_backup.FrameArray(i,1) + AOA_backup.FrameArray(i,2) ) / 2 );
        
        if Lt > 1 % (1.2)
            DISPLAY.step = i;
        end
        
        %%% For each Display to be plot
        for j = 1:nQname
            
            if isfield(AOA_backup,Qname{j}) % 2.3
                
                % Get the Qcolor, Qunits, Qsr, etc. corresponding to the Display
                [Pname,idx] = GetPname(Qname{j}); % now gets Pname (2.3)
                Qsr       = eval(['sr_' Pname '{idx}']);
                Qscalebar = eval(['srbar_' Pname '(idx)']);
                Qunits    = eval(['allUnits_' Pname '{idx}']);
                Qcolor    = eval(['allColors_' Pname '{idx}']);
                QKillTr   = eval(['killtrace_' Pname '(idx)']);
                
                % tag indicating mean Tr = 0 for naming files
                KillTrTag = '';
                if QKillTr
                    KillTrTag = '_Tr=0';
                end

                % Load significance map
                if PLOT.significance
                    temp_sign_map = eval(['AOA_backup.' Qname{j} '_Smap']);
                    DISPLAY.significant_map = temp_sign_map(:,:,:,i);
                end
                
                % Plot image
                PlotField(Qname{j}, QKillTr, Qcolor , Qunits , Qsr, Qscalebar, AOA_backup, Qimage, DISPLAY);
                                
                % Plot print (removed "Frames_" prefix and added fullTag in filenames)(2.6)
                if PLOT.print
                    thisFilename = ['AOA_' Qname{j} '_' timeStart 'to' timeEnd KillTrTag fullTag '_sr=' num2str(Qsr(1)) imageExtension]; % 1.7, added fullTag (2.5), use timeStart/End (2.7)
                    if ~exist([OutputPathName filesep 'AOA_' QplotType fullTag],'dir')
                        mkdir([OutputPathName filesep 'AOA_' QplotType fullTag]); 
                    end
                    if strcmp(imageExtension,'.svg')
                        plot2svg([OutputPathName filesep 'AOA_' QplotType fullTag filesep thisFilename],figure(1),'png');
                    else
                        print(printFormat, printResolution, [OutputPathName filesep 'AOA_' QplotType fullTag filesep thisFilename]);
                    end
                end
                close
                
            else % 2.3
                disp(['AOA WARNING: quantity "' Qname{j} '" was not found in "AOA_backup"!']);
                disp('Please make sure that it is the name assigned to this quantity AND that it was averaged over time and stored in each animal AOT backups.')
            end        
        end % end for each quantities
    end % end for each time point
end % end if plot


%% History %%

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




