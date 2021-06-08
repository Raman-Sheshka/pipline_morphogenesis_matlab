% AverageOverAnimals (AOA)
%
% Script containing process to average data from multiple AOT runs
%
% Version 2.6
% Stephane Rigaud
% Anais Bailles
% Boris Guirao


%% Building of global grid %%

% Merge all the GRID information loaded from the different annimals in mean_GRID_specs
% Needed for determining the grid on which to mean and the plot are made.

disp('Build global grid ...')

% Load grids to define crop and translation
clear GRID;
Ox = zeros(size(allAnimals)); % use to be size but it bug now ... ???
Oy = zeros(size(allAnimals)); % use to be size but it bug now ... ???
lx = zeros(size(allAnimals)); % use to be size but it bug now ... ???
ly = zeros(size(allAnimals)); % use to be size but it bug now ... ???

for n = 1:nAllAnimals % for all the animal
    grid = mean_GRID_specs{n};
    
    % find the coordinates (integers) in box number of the barycenter
    ULC = cell2mat3D(grid.ULCs);
    [useless1 Ox_a] = ind2sub(grid.size, find( ULC(:,:,1) == grid.xywh(1)));
    [Oy_a useless2] = ind2sub(grid.size, find( ULC(:,:,2) == grid.xywh(2)));
    Ox(n) = unique(Ox_a); % number of boxes (along x) starting from left side of grid
    Oy(n) = unique(Oy_a); % number of boxes (along y) starting from top side of grid
    
    % find number of compartments along x and y
    lx(n) = grid.size(2);
    ly(n) = grid.size(1);
    
    % if movie_type(n) == 2 % if this is left part
    % Oy(n) = ly(n)-(Oy(n)-1)+1; % coordinates of the origine once the grid has been flipped with Y(Rflipped)=Ly-Y(Runflipped)+1
    % end
end

% Find total coverture grid size
Ltotx = max(Ox)+max(lx-Ox); % lx-Ox = number of boxes at the right of box Ox
Ltoty = max(Oy)+max(ly-Oy); % ly-Oy = number of boxes below box Oy

% Spot grid location in global common grid for each animal
startLocXwt = zeros(nAvgAnimals,1);
startLocYwt = zeros(nAvgAnimals,1);
endLocXwt = zeros(nAvgAnimals,1);
endLocYwt = zeros(nAvgAnimals,1);
for n = 1:nAvgAnimals
    [~,idx] = ismember(avgAnimals{n}, allAnimals);
    startLocXwt(n) = max(Ox) - Ox(idx) + 1; % +1 stems from the location being called from 1 to n (cf the piquets and barrieres pb)
    endLocXwt(n) = max(Ox) - Ox(idx) + lx(idx);
    startLocYwt(n) = max(Oy) - Oy(idx) + 1;
    endLocYwt(n) = max(Oy) - Oy(idx) + ly(idx);
end

% determine the image size that can contain the union grid, based on the box size
% and the new origin position in the image, for display purposes
sizeImageX = ceil(PLOT.boxSize(1) * Ltotx * (1 - gridOverlap))*1.05;     % added "ceil" (2.3)
sizeImageY = ceil(PLOT.boxSize(2) * Ltoty * (1 - gridOverlap)); % added "ceil" and leeway of 5% (2.3)
newOx = max(Ox) * PLOT.boxSize(1) * (1 - gridOverlap) - PLOT.boxSize(1)/2;
newOy = max(Oy) * PLOT.boxSize(2) * (1 - gridOverlap) - PLOT.boxSize(2)/2;
Qimage = ones([sizeImageY , sizeImageX]);

% Build the global common grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRID = GridMaker([sizeImageY sizeImageX], PLOT.boxSize, [newOx newOy], [] , gridColor, gridLineWidth, gridOverlap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
macrocaetes = DISPLAY.macrocaetes;
meanYmid = newOy - meanDeltaYmid; % 2.5

%%% Defining FULL animal quantities (mod 2.6):
% Redefining newOy, sizeImageY, GRID and macrocaetes if making full animal out of halves (2.5)
extraImageY = sizeImageY - meanYmid;        % will only flip the part beyond midline
newOyFull = newOy + extraImageY;
meanYmidFull = newOyFull - meanDeltaYmid;           % meanYmid update
sizeImageYFull = 2*sizeImageY;
QimageFull = ones([sizeImageYFull , sizeImageX]);
gridSizeFull = GRID.size;                       % gets original grid size that will be kept
% since using gridSize to keep previous size, need to determine grid ew ULC in extended image
gridULCx = PLOT.boxSize(1)/2;
gridULCy = PLOT.boxSize(2)/2 + extraImageY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIDfull = GridMaker([sizeImageYFull sizeImageX], PLOT.boxSize, [gridULCx gridULCy], gridSizeFull , gridColor, gridLineWidth, gridOverlap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIDfull.xywh = [newOx newOyFull GRID.xywh(3) GRID.xywh(4)];
macrocaetesFull(1,:) = macrocaetes(1,:) + newOx;
macrocaetesFull(2,:) = macrocaetes(2,:) + newOyFull;
% filling coordinates for 8 flipped macrocaetes:
macrocaetesFull(1,9:16) = macrocaetes(1,1:8);
macrocaetesFull(2,9:16) = 2*meanYmidFull - macrocaetes(2,1:8); % (yMacro + yMacroFlipped)/2 = meanYmid => yMacroFlipped = ...

%%% Regular (not full) case:
macrocaetes(1,:) = macrocaetes(1,:) + newOx;
macrocaetes(2,:) = macrocaetes(2,:) + newOy;

fullTag = '';
if PLOT.makeItFull
    fullTag = '_full';
    DISPLAY.macrocaetes = macrocaetesFull; % update
    DISPLAY.yMid = meanYmidFull; % 2.5
else
    DISPLAY.macrocaetes = macrocaetes; % update
    DISPLAY.yMid = meanYmid; % 2.5
end




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

% AOA backup initialisation (mod 2.6)
% We try to mimic the structure of single animal backup
if PLOT.makeItFull
    AOA_backup = GRIDfull;
    AOA_stdMAP = GRIDfull;
    QimagePlot = QimageFull;
else
    AOA_backup = GRID;
    AOA_stdMAP = GRID;
    QimagePlot = Qimage;
end

% find the number of time point
[Lt] = size( eval(['GRID_AOT_' avgAnimals{1} '.' list{1}]) , 4 ); % 2.3

%%% Stackin up quantities
%----------------------------------------------------------------------------------------
% for each quantities, including AreaRatios and Rcons, we create a 5D-matrix (X,Y,Z,Time,Animal)
% wherw X and Y are defined by the new GRID, Z the dimension of the quantity values (1 scalar, 2 vector, or 4 tensor), Time the number of time point and Animal the number of animal
fprintf('Stacking up data ...')
for m = 1:nList
    [~,~,Lz,~] = size( eval(['GRID_AOT_' avgAnimals{1} '.' list{m}]) ); % get the dimensions of the quantity (2.3)
    eval(['total_' list{m} '= NaN(Ltoty, Ltotx, Lz, Lt, size(avgAnimals,1));']); % create a 5D matrix
    
    % Specific case for AreaRatios (2.3)
    if strcmp(list{m},'AreaRatios')
        total_AreaRatios = zeros(Ltoty, Ltotx, Lz, Lt, size(avgAnimals,1)); % NB: RConds already included in AreaRatios for TA
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
if exist('macrocaetes','var')
    AOA_backup.Macrocaetes = DISPLAY.macrocaetes;
    AOA_stdMAP.Macrocaetes = DISPLAY.macrocaetes; % (v2) in progress add-on for biological variability estimation
end

%%% Save AOA mean backup
if ~exist(OutputPathName,'dir'); mkdir(OutputPathName); end;
save([OutputPathName filesep 'AOA_stdMAP'], 'AOA_stdMAP');
if ~exist([OutputPathName filesep 'AOA_backup.mat'],'file')
    disp('Saving new backup ...')
    save([OutputPathName filesep 'AOA_backup'], '-struct', 'AOA_backup');
else
    disp('Appending backup ...')
    Previous = load([OutputPathName filesep 'AOA_backup']);
    AOA_backup = catstruct(Previous,AOA_backup);
    save([OutputPathName filesep 'AOA_backup'], '-struct', 'AOA_backup', '-append');
end


%% Plot %%

if ~isempty(Qname) && PLOT.plot
    
    disp('Starting plot ...')
    
    AOA_backup = load([OutputPathName filesep 'AOA_backup']);
%     AOA_backup = load([OutputPathName filesep 'AOA_backup_MOD']);
    % NB: "MOD" backups were made to symmetrize tensors ONLY in the row right next to the midline, 
    % to mimick the calculation done in BIGwt2 were a big box (similar to Jesus' box) including
    % both sides of the midline can be drawn.

    %%% Initialisation of the plot
    DISPLAY.minimalInfoDisplay = minimalInfoDisplay;
    DISPLAY.scaleBarWidth = scaleBarWidth;
    DISPLAY.gridDisplay = gridDisplay;
    DISPLAY.lineWidth = lineWidth;
    DISPLAY.pointSize = pointSize;
    DISPLAY.Animal = Animal;
    DISPLAY.EVstyles = EVstyles;
    DISPLAY.signOpacities = signOpacities;
    DISPLAY.fontSizeInfo = fontSize;
    DISPLAY.imageFading = imageFading;
    DISPLAY.errorPsMin = 1;
    DISPLAY.errorDnPsMin = 1;
    DISPLAY.errorFontSize = 5;
    DISPLAY.plotType = QplotType;
    DISPLAY.fadeColor = custom_white; % dont know why but plotfield need it for VM
    DISPLAY.displayOrigin = true;
    DISPLAY.drawMidline = PLOT.drawMidline; % 2.5
    
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
                
                if Lt > 1 % (1.2)
                    DISPLAY.step = i;
                end
              
                % Load significance map
                if PLOT.significance
                    temp_sign_map = eval(['AOA_backup.' Qname{j} '_Smap']);
                    DISPLAY.significant_map = temp_sign_map(:,:,:,i);
                end
                
                % Plot image
                PlotField(Qname{j}, QKillTr, Qcolor , Qunits , Qsr, Qscalebar, AOA_backup, QimagePlot, DISPLAY); % use of QimagePlot (2.6)
                                
                % Plot print (removed "Frames_" prefix to folder names  2.5)
                if PLOT.print
                    thisFilename = ['AOA_' Qname{j} '_' AOA_backup.TimeArray{i,1} 'to' AOA_backup.TimeArray{i,2} KillTrTag fullTag '_sr=' num2str(Qsr(1)) imageExtension]; % 1.7, added fullTag (2.5)
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

% 20/01/2017: 2.6 (Boris)

% 18/01/2017: 2.5 (Boris)

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




