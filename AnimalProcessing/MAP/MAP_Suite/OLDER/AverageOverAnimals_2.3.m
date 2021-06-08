% AverageOverAnimals (AOA)
%
% Script containing process to average data from multiple AOT runs
%
% Version 2.3
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
    
    % find the coordinates in box number of the barycenter
    ULC = cell2mat3D(grid.ULCs);
    [useless1 Ox_a] = ind2sub(grid.size, find( ULC(:,:,1) == grid.xywh(1)));
    [Oy_a useless2] = ind2sub(grid.size, find( ULC(:,:,2) == grid.xywh(2)));
    Ox(n) = unique(Ox_a);
    Oy(n) = unique(Oy_a);
    
    % find length
    lx(n) = grid.size(2);
    ly(n) = grid.size(1);
    
    % if movie_type(n) == 2 % if this is left part
    % Oy(n) = ly(n)-(Oy(n)-1)+1; % coordinates of the origine once the grid has been flipped with Y(Rflipped)=Ly-Y(Runflipped)+1
    % end
end

% Find total coverture grid size
Ltotx = max(Ox)+max(lx-Ox);
Ltoty = max(Oy)+max(ly-Oy);

% Spot grid location in total coverture grid for each animal
locx_wt = zeros(size(avgAnimals));
locy_wt = zeros(size(avgAnimals));
endlocx_wt = zeros(size(avgAnimals));
endlocy_wt = zeros(size(avgAnimals));
for n = 1:nAvgAnimals
    [~,idx] = ismember(avgAnimals{n}, allAnimals);
    locx_wt(n) = max(Ox) - Ox(idx) + 1; % +1 stems from the location being called from 1 to n (cf the piquets and barri?res pb)
    endlocx_wt(n) = max(Ox) - Ox(idx) + lx(idx);
    locy_wt(n) = max(Oy) - Oy(idx) + 1;
    endlocy_wt(n) = max(Oy) - Oy(idx) + ly(idx);
end

% determine the image size that can contain the union grid, based on the box size
% and the new origin position in the image, for display purposes
sizeImageX = ceil(PLOT.boxSize(1) * Ltotx * (1 - gridOverlap));     % added "ceil" (2.3)
sizeImageY = ceil(PLOT.boxSize(2) * Ltoty * (1 - gridOverlap))*1.1; % added "ceil" and leeway of 10% (2.3)
newOx = max(Ox) * PLOT.boxSize(1) * (1 - gridOverlap) - PLOT.boxSize(1)/2;
newOy = max(Oy) * PLOT.boxSize(2) * (1 - gridOverlap) - PLOT.boxSize(2)/2;
Qimage = ones([sizeImageY , sizeImageX]);

% %-----------
% % determine the image size that can contain the union grid, based on the box size
% % and the new origin position in the image, for display purposes
% % NB: specific case if grid overlap is 0
% size_image_x = PLOT.boxSize(1) * (Ltotx * gridOverlap + 1);
% size_image_y = PLOT.boxSize(2) * (Ltoty * gridOverlap + 1);
% newOx = max(Ox) * PLOT.boxSize(1) * gridOverlap + 1;
% newOy = max(Oy) * PLOT.boxSize(2) * gridOverlap + 1;
% if gridOverlap == 0
%     size_image_x = PLOT.boxSize(1) * (Ltotx + 1);
%     size_image_y = PLOT.boxSize(2) * (Ltoty + 1);
%     newOx = max(Ox) * PLOT.boxSize(1) + 1;
%     newOy = max(Oy) * PLOT.boxSize(2) + 1;
% end
% Qimage = ones([size_image_y , size_image_x]);
% %-------------------

% build the union grid
GRID = GridMaker([sizeImageY sizeImageX], PLOT.boxSize, [newOx newOy], [] , gridColor, gridLineWidth, gridOverlap);
DISPLAY.macrocaetes(1,:) = DISPLAY.macrocaetes(1,:) + newOx;
DISPLAY.macrocaetes(2,:) = DISPLAY.macrocaetes(2,:) + newOy;


%% Initialization %%

% set up and test all variable for the averaging process
% NB: specific case if averaging AOS process for managin multichannel

fprintf('Uniformizing data format ...')

% signal_prefixes = {'CD_' 'mCD_'};  % possible expected prefix in multichannel process
% eval(['list = ' PnameAOA '_list;']);  % load quantities list associated to the process

eval( ['temporary = GRID_AOT_' avgAnimals{1} ';' ] ); % 2.3
% eval( ['temporary = GRID_' PnameAOA '_' avgAnimals{1} ';' ] );
list = fieldnames(temporary);
for m = 2:nAvgAnimals
    eval( ['temporary = GRID_AOT_' avgAnimals{m} ';' ] ); % 2.3
%     eval( ['temporary = GRID_' PnameAOA '_' avgAnimals{m} ';' ] );
    list = intersect(list,fieldnames(temporary));
end
list = setdiff(list,rejectList);
nList = length(list); % 2.3

% dimension test
eval( ['nbFrame = size(GRID_AOT_' avgAnimals{1} '.FrameArray,1 );' ] );% 2.3
% eval( ['nbFrame = size(GRID_' PnameAOA '_' avgAnimals{1} '.FrameArray,1 );' ] );
testDim = 4;
if nbFrame == 1, testDim = 3; end

%%% Turning all 3D matrix (X,Y,Time) into 4D matrix (X,Y,1,Time)
for n = 1:nList
    for m = 1:nAvgAnimals
        eval( ['temporary = GRID_AOT_' avgAnimals{m} '.' list{n} ';' ] ); % 2.3
%         eval( ['temporary = GRID_' PnameAOA '_' avgAnimals{m} '.' list{n} ';' ] );
        if size(temporary,testDim) == 1
            new = reshape( temporary, size(temporary,1), size(temporary,2) ,1, size(temporary,3) );
            eval( ['GRID_AOT_' avgAnimals{m} '.' list{n} ' = new ;' ] ); % 2.3
%             eval( ['GRID_' PnameAOA '_' avgAnimals{m} '.' list{n} ' = new ;' ] );
        end
        % OLD
%         signal_TF = strncmp(list{n}, signal_prefixes, 2); % returns logical with 1 where 'CD_' or 'mCD_' were found
%         if any(signal_TF)
%             eval( ['temporary = GRID_' PnameAOA '_' avgAnimals{m} '.' list{n} AOS_dif_name{m} ';' ] );
%             if size(temporary,testDim) == 1
%                 new = reshape( temporary, size(temporary,1), size(temporary,2) ,1, size(temporary,3) );
%                 eval( ['GRID_' PnameAOA '_' avgAnimals{m} '.' list{n} AOS_dif_name{m} ' = new ;' ] );
%             end
%         else
%             eval( ['temporary = GRID_' PnameAOA '_' avgAnimals{m} '.' list{n} ';' ] );
%             if size(temporary,testDim) == 1
%                 new = reshape( temporary, size(temporary,1), size(temporary,2) ,1, size(temporary,3) );
%                 eval( ['GRID_' PnameAOA '_' avgAnimals{m} '.' list{n} ' = new ;' ] );
%             end
%         end
    end % end for each animal
end % end for each quantities
fprintf('Done\n')

disp('FYI: Loaded backup to be averaged ...')
for m = 1:nAvgAnimals
    disp(backupPathList{m}); % display the backup being loaded
    %disp(['(' num2str(locx_wt(m)) ',' num2str(locy_wt(m)) ') ('  num2str(endlocx_wt(m)) ',' num2str(endlocy_wt(m)) ')']);
end


%% Calculating averages over animals %%

% AOA backup initialisation
% We try to mimic the structure of single animal backup
AOA_backup = GRID;
AOA_stdMAP = GRID;

% find the number of time point
[Lt] = size( eval(['GRID_AOT_' avgAnimals{1} '.' list{1}]) , 4 ); % 2.3
% [Lt] = size( eval(['GRID_' PnameAOA '_' avgAnimals{1} '.' list{1}]) , 4 );

%%% Stackin up quantities
%----------------------------------------------------------------------------------------
% for each quantities, including AreaRatios and Rcons, we create a 5D-matrix (X,Y,Z,Time,Animal)
% wherw X and Y are defined by the new GRID, Z the dimension of the quantity values (1 scalar, 2 vector, or 4 tensor), Time the number of time point and Animal the number of animal
fprintf('Stacking up data ...')
for m = 1:nList
    [~,~,Lz,~] = size( eval(['GRID_AOT_' avgAnimals{1} '.' list{m}]) ); % get the dimensions of the quantity (2.3)
%     [~,~,Lz,~] = size( eval(['GRID_' PnameAOA '_' avgAnimals{1} '.' list{m}]) ); % get the dimensions of the quantity
    % OLD
%     signal_TF = strncmp(list{m}, signal_prefixes, 2); % returns logical with 1 where 'CD_' or 'mCD_' were found
%     if any(signal_TF)
%         [~,~,Lz,~] = size( eval(['GRID_' PnameAOA '_' avgAnimals{1} '.' list{m} AOS_dif_name{1}]) ); % (AOS compatibility)
%     else
%         [~,~,Lz,~] = size( eval(['GRID_' PnameAOA '_' avgAnimals{1} '.' list{m}]) ); % get the dimensions of the quantity
%     end

    eval(['total_' list{m} '= NaN(Ltoty, Ltotx, Lz, Lt, size(avgAnimals,1));']); % create a 5D matrix
    
    % Specific case for AreaRatios (2.3)
    if strcmp(list{m},'AreaRatios')
        total_AreaRatios = zeros(Ltoty, Ltotx, Lz, Lt, size(avgAnimals,1)); % NB: RConds already included in AreaRatios for TA
    end
    
%     % Specific cases for AreaRatios & RConds
%     if strncmp(list{m},'AreaRatios',10)
%         total_AreaRatios = zeros(Ltoty, Ltotx, Lz, Lt, size(avgAnimals,1)); 
%     elseif strcmp(list{m},'RConds')
%         total_RConds = zeros(Ltoty, Ltotx, Lz, Lt, size(avgAnimals,1));
%     end
    
    for a = 1:nAvgAnimals
        eval(['total_' list{m} '(locy_wt(a):endlocy_wt(a), locx_wt(a):endlocx_wt(a),:,:,a) = GRID_AOT_' avgAnimals{a} '.' list{m} ';']); % 2.3
        % OLD
%         if any(signal_TF)
%             eval(['total_' list{m} '(locy_wt(a):endlocy_wt(a), locx_wt(a):endlocx_wt(a),:,:,a) = GRID_' PnameAOA '_' avgAnimals{a} '.' list{m} AOS_dif_name{a} ';']);
%         else
%             eval(['total_' list{m} '(locy_wt(a):endlocy_wt(a), locx_wt(a):endlocx_wt(a),:,:,a) = GRID_' PnameAOA '_' avgAnimals{a} '.' list{m} ';']);
%         end
    end % end for each animal
end % end for each quantities
fprintf('Done\n')
%----------------------------------------------------------------------------------------

%%% Calculating weights
%----------------------------------------------------------------------------------------
fprintf('Calculating weights ...')
total_weight = total_AreaRatios .^ 2;
% if exist('total_RConds','var')
%     total_weight = (total_AreaRatios .* total_RConds) .^ 2;
%     total_weight(isnan(total_weight)) = 0;
% end
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
        
%     if strncmp(list{m},'AreaRatios',10)
%         disp(['    Processing ' list{m} '...']);
%         eval(['AOA_backup.AreaRatios_' PnameAOA '= mean( total_AreaRatios, 5 );'])
%         eval(['AOA_stdMAP.AreaRatios_' PnameAOA '= mean( total_AreaRatios, 5 );'])
%         
%     elseif strcmp(list{m},'RConds')
%         disp(['    Processing ' list{m} '...']);
%         AOA_backup.RConds = nanmean( total_RConds, 5 );  % nanmean because, RConds should not impact weight when not calculated
%         AOA_stdMAP.RConds = nanmean( total_RConds, 5 );
%         tmpIndex = AOA_backup.RConds == 0;               % nanmean of a full NaN array return 0, we switch back 0s to NaN
%         AOA_backup.RConds(tmpIndex) = NaN;
%         AOA_stdMAP.RConds(tmpIndex) = NaN;
        
    else
        disp(['    Processing ' list{m} '...']);
        % repeat weight matrix to fit the quantity structure (scalar, vector, or tensor)
        [Tz] = eval(['size( total_' list{m} ',3);']);
        R_total_weight = repmat(total_weight, [1 1 Tz 1 1]);
        % calculate the weighted average of the quantity
        temp_mean = eval(['total_' list{m} ' .* R_total_weight;']);
        temp_mean = nansum( temp_mean, 5) ./ sum( R_total_weight, 5 );
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
% if ~exist(OutputPathName,'dir'); mkdir(OutputPathName); end;
% save([OutputPathName filesep 'AOA_' PnameAOA '_backup'], ['AOA_' PnameAOA '_backup']);


%% Plot %%

if ~isempty(Qname) && PLOT.plot
    
    disp('Starting plot ...')
    
    AOA_backup = load([OutputPathName filesep 'AOA_backup']);
%     AOA_backup = load([OutputPathName filesep 'AOA_backup_MOD']);
    % NB: "MOD" backups were made to symmetrize tensors ONLY in the row right next to the midline, 
    % to mimick the calculation done in BIGwt2 were a big box (similar to Jesus' box) including
    % both sides of the midline can be drawn.

    %%% Initialisation of the plot
%     signal_prefixes = {'CD_' 'mCD_'};
    
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
    
    %%% For each time point
    for i = 1:Lt
        
        % Set the time and n variable (mod 2.2)
        timeStart = Time_str2dec(AOA_backup.TimeArray{i,1}) + timeShift;                % turns time string into decimal and adds timeShift
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
                
                % load the coresponding AreaRatios
%               eval(['AOA_backup.AreaRatios = AOA_backup.AreaRatios_' PnameAOA ';']);
                
                % Load significance map
                if PLOT.significance
                    temp_sign_map = eval(['AOA_backup.' Qname{j} '_Smap']);
                    DISPLAY.significant_map = temp_sign_map(:,:,:,i);
                end
                
                % Plot image
                PlotField(Qname{j}, QKillTr, Qcolor , Qunits , Qsr, Qscalebar, AOA_backup, Qimage, DISPLAY);
                
                % Plot print
                if PLOT.print
                    this_filename = ['AOA_' Qname{j} '_' AOA_backup.TimeArray{i,1} 'to' AOA_backup.TimeArray{i,2} KillTrTag '_sr=' num2str(Qsr(1)) imageExtension]; % 1.7
                    if ~exist([OutputPathName filesep 'Frames_AOA_' QplotType],'dir'); mkdir([OutputPathName filesep 'Frames_AOA_' QplotType]); end;
                    if strcmp(imageExtension,'.svg')
                        plot2svg([OutputPathName filesep 'Frames_AOA_' QplotType filesep this_filename],figure(1),'png');
                    else
                        print(printFormat, printResolution, [OutputPathName filesep 'Frames_AOA_' QplotType filesep this_filename]);
                    end
                end
                close
                
            else % 2.3
                disp(['AOA WARNING: quantity "' Qname{j} '" was not found in "AOA_backup"!']);
                disp('Please make sure that it is the name assigned to this quantity AND that it was averaged over time and stored in each animal AOT backups.')
            end
            
            %% OLD
%             % 'CD' and 'mCD' cases management if AOA over 'AOS'
%             signal_TF = strncmp(Qname{j}, signal_prefixes, 2); % returns logical with 1 where 'CD_' or 'mCD_' were found
%             if any(signal_TF) % (AOS compatibility)
%                 for r = 1:length(AOS_sub_name)
%                     % Load significance map
%                     if PLOT.significance
%                         temp_sign_map = eval(['AOA_backup.' Qname{j} '_' AOS_sub_name{r} '_Smap']);
%                         DISPLAY.significant_map = temp_sign_map(:,:,:,i);
%                     end
%                     
%                     % Plot image
%                     PlotField([Qname{j} '_' AOS_sub_name{r}], QKillTr, Qcolor, Qunits , Qsr, Qscalebar, AOA_backup, Qimage, DISPLAY);
%                     
%                     % Plot print
%                     if PLOT.print
%                         this_filename = ['AOA_' Qname{j} '_' AOS_sub_name{r} '_' AOA_backup.TimeArray{i,1} 'to' AOA_backup.TimeArray{i,2} KillTrTag '_sr=' num2str(Qsr(1))  image_extension_output]; % 1.7
%                         if ~exist([OutputPathName filesep 'Frames_AOA_' QplotType],'dir'); mkdir([OutputPathName filesep 'Frames_AOA_' QplotType]); end;
%                         if strcmp(image_extension_output,'.svg')
%                             plot2svg([OutputPathName filesep 'Frames_AOA_' QplotType filesep this_filename],figure(1),'png');
%                         else
%                             print(print_format, print_resolution, [OutputPathName filesep 'Frames_AOA_' QplotType filesep this_filename]);
%                         end
%                     end
%                     close
%                 end % end for each channels
%             else
%                 % Load significance map
%                 if PLOT.significance
%                     temp_sign_map = eval(['AOA_backup.' Qname{j} '_Smap']);
%                     DISPLAY.significant_map = temp_sign_map(:,:,:,i);
%                 end
%                 
%                 % Plot image
%                 PlotField(Qname{j}, QKillTr, Qcolor , Qunits , Qsr, Qscalebar, AOA_backup, Qimage, DISPLAY);
%                 
%                 % Plot print
%                 if PLOT.print
%                     this_filename = ['AOA_' Qname{j} '_' AOA_backup.TimeArray{i,1} 'to' AOA_backup.TimeArray{i,2} KillTrTag '_sr=' num2str(Qsr(1)) imageExtension]; % 1.7
%                     if ~exist([OutputPathName filesep 'Frames_AOA_' QplotType],'dir'); mkdir([OutputPathName filesep 'Frames_AOA_' QplotType]); end;
%                     if strcmp(imageExtension,'.svg')
%                         plot2svg([OutputPathName filesep 'Frames_AOA_' QplotType filesep this_filename],figure(1),'png');
%                     else
%                         print(printFormat, printResolution, [OutputPathName filesep 'Frames_AOA_' QplotType filesep this_filename]);
%                     end
%                 end
%                 close
%             end 
            
        end % end for each quantities
    end % end for each time point
    
end % end if plot


%% History %%

% 12/01/2017: 2.3 (Boris)
% - ONLY using AOT backups now!! (=> removed "PnameAOA")

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




