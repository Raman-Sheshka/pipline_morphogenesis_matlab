% Advanced_Image_Analysis (AIA)
%
% Version 2
% Stéphane Rigaud
%
% Script containing all necessary info required to run AOA ... and is aimed to extract most of thos information respectively
% from the different AIA_Info of each animals.

% definitions and conventions
% Quantities : EG, ED, EA, I, V, S, ... ; shortened into Q
% Processes  : AOS, TA, SM, ... ; shortened into P
% Operations : AOA, DOA, ... ; shortened into O
% 
% With exception of AOA, all operations are first intended to take operation 
% backup in input (specified by the Oname variable)
% If Oname is left empty, the operation will expect a single animal


%% DO NOT MODIFY
% Sanity cleaning of Matlab
clear all; close all; clc; 
% Flag preventing AIA_parameters to run other programs
AIA_call = 1; 


%% Operation Name (Oname) Selection
% Boolean in order to select which Operation will be run ; 
% NB: /!\ only run one at a time

AOA = 1;               % Average Over Animal
DOA = 0;               % Difference Over Animal
COQ = 0;               % Contraction Over Quantities (in a same animal)
COA = 0;               % Contraction Over Animal
AOZ = 0;               % Analysis Over Zones (temporal plot)

TC  = 0;               % Tensor Corelator (in progress)

% Global output path. All the calculated output will be saved in a file
% structure located at this specified path.
PathName = 'D:\partage\AOA_Outputs_test';
AIAFolderName = 'C:\Users\Stephane\Documents\MATLAB\final\v2_2016-03-15'; % path to the folder containing the AIA_info_X
 


%% Common Parameters
% Similar parameters as in AIA_Parameters.
CustomColors;       % defines usual set of colors
AllQsColorsUnits;   % associate quantities with specific colors

% display parameters to be used during plot
minimalInfoDisplay = false;
minAEV = 0.0001;
fontSize = 12;
EVstyles = {'-' '--'};       % ONLY relevant for "merged" display type: styles to display ellipse axes representing tensor eigenvalues (default {'-' ':'})
pointSize = 2;
signOpacities = [0.7 0.3];   % ONLY relevant for "split+/-" and "circle" display types: specifies opacity of positive(white) and negative(black) disks, respectively.
lineWidth = 1.5;             % for circle, bars and ellipses (1.5 ok with BIG movies)
gridDisplay = false;         % Lagrangian grid ALWAYS displayed (6.0)
gridColor = black;
gridLineWidth = 0.5;        % only matters for Egrid, Lgrid thickness specified in LGridPlotter
imageFading = 0.6;
scaleBarWidth = 1;

% plot process parameters
PLOT.print = 1;                  % print the plot
PLOT.plot = 1;                   % plot the selected quantities (Qname)
PLOT.significance = 1;           % plot significance
PLOT.SignOpacityMap = [0.5 0.5]; % significance opacity [iso dev]
PLOT.extension = 'png';          % plot image type (png | svg)
PLOT.resolution = 300;           % plot resolution
PLOT.boxSize = [128 128];        % box size of the average grid
PLOT.macrocaetes = 1;            % plot macrocaetes position
PLOT.ARtreshold = 0;             % AreaRatios treshold (in progress)

%% Scaleratio, Scalebar, and Trace values
% Scale values according to the information to be ploted
sr_AOS        = { 8e2 ; 50 ;  2 ; 100 ; [4 30] };   % sets ratio setting size of ellipses or bars in tensor representation for M,I.
srbar_AOS     = [ 0.1 ;  2 ; 50 ;  1  ;    10  ];   % scale bar lengths for each contribution
killtrace_AOS = [  0  ;  0 ;  0 ;  1  ;     1  ];   % Will set average compartment trace to 0 in the plots (mean isotropic part = 0). Choose this when tensors are known up to an additive constant
% contributions   Rho    I    M    V       CD

sr_SM        = { 1500 ; 1000 ; 1000 ; 1000 };
srbar_SM     = [  0.1 ; 0.05 ; 0.05 ; 0.05 ];
killtrace_SM = [    1 ;    1 ;    1 ;    1 ];
% contributions     S     SP     ST      P

sr_TA        = {  4e3 ;  4e3 ;  4e3 ;  8e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 };  % average 14h
srbar_TA     = [ 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ];
killtrace_TA = [    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ];
% contributions     G      S      R     Ds      D      A      N      F      J     Jb     DM      U     G*    PSI   PhiU

animalTimeWidth   = 13;
animalTimeOverlap = 0.5;
multiTimeStart    = '14h55'; 
multiTimeStop     = '27h55';

%% AOA - Average over Animals
if AOA  
    % This process calculate the average of N individual (N can be 1) together. 
    % The spatial and temporal average should be the same between the animal 
    % in order to average them together (see AOT)
    
    % Specify the list of all animals name that will be use in the generation of the average grid
    % NB: The list should contain all animal that will be compared and not only the animal to be ploted
    % by default: provide all the possible animals
    All_Animals = { 'BIG_1' ; 'BIG_5' ; 'BIG_6' ; 'BIGwt2r' ; 'BIGwt2l' ; ...
                    'TRBL4' ; 'TRBL7' ; 'TRBL9' ; 'TRBL8r'  ; 'TRBL8l'  }; % animal name list
    
    genericName = 'single';
    
    %%% Values
    % Pname = | SM ; VMM ; AOS ; TA | ... 
    % Qname = | EA ; EG ; ED ; mI ; etc. | ... 
    % QPlot_type = | dev+/- ;  split+/- ; merge |
    
    Pname = 'TA';                                             % Process name
    Qname = {'EA';'ER';'ES';'ED';'EDs';'EJ';'EN';'EF';'EG'};  % Quantities name
    % QplotType = 'dev+';                                       % Plot Type
    % QplotType = 'circle';
    QplotType = 'split+';
    
    % Animals = { 'BIG_1' };
    Animals = { 'BIG_1' ; 'BIG_5' ; 'BIG_6' ; 'BIGwt2r' ; 'BIGwt2l'}; % list of animals to be processed
    % Animals = { 'TRBL4' ; 'TRBL7' ; 'TRBL9' ; 'TRBL8r'  ; 'TRBL8l' };
    macrocaete_folder = 'C:\Users\Stephane\Desktop\10-12-2014\output_trbl\Landmarks_output_2_R';
end


%% DOA - Difference over Animals
if DOA
    % This process calculate the difference between two individual.
    % In this case, an individual is the output of a Process (AOA, DOA, COA, other)
    % In the current state, both individual are coming from the same Process (AOA - AOA, DOA - DOA, but no AOA - DOA)
    % As this is a simple A - B operation, both individual must have the same spatial and temporal size, with same overlap
    
    Oname = 'AOA'; % Operation name, single animal process if empty
    
    % Animal 1 minus Animal 2
    % The name of the output animal will be A1-A2
    Animal  = { 'BIG_1' ; 'BIG_2'};
    DOAname = [Animal{1} '-' Animal{2}];
    
    %%% Values
    % Pname = | SM ; VMM ; AOS ; TA | ... 
    % Qname = | EA ; EG ; ED ; mI ; etc. | ... 
    % QPlot_type = | dev+/- ;  split+/- ; merge |
    
    Pname = 'TA';                                             % Process name
    Qname = {'EA';'ER';'ES';'ED';'EDs';'EJ';'EN';'EF';'EG'};  % Quantities name
    QplotType = 'dev+';                                       % Plot Type
    % QplotType = 'circle';
    
    PLOT.Qparameters = 1;
    % if Qparameters is off, use the following parameters for the plot
    Qcolor = [0.3 0.5 0.5];
    Qunits = 'A.U.';
    Qsr = 4e3;
    Qsrbar = 2e-2;
    QKillTr = 0;  
end


%% COQ - Calcul over Quantities
if COQ
    
    Animal = 'big';
    
    Oname = 'AOA';
    
    Qname1 = {'ED'};
    Qname2 = 'EG';

    % if Qparameters is off, use the following Qparameters
    Qcolor = [0.3 0.5 0.3];
    Qunit = 'A.U.';
    Qsr = 4e3;
    Qsrbar = 2e-2;
    QKillTr = 0;
    
    QplotType = 'dev+';
end


%% COA - Contraction over Animals
if COA  %%% TODO: test on single animal ex: BIGwt2, lunch from AIA_paramters
        
    tagPlot = 'd';  % d i do
    
    %%% Process parameters
    unitaire = 0;
    wtContraction = 0;
    
    globalTime = 17.25; globalOverlap = 0;
    
    Oname = 'DOA';
    Animal = 'big-trbl';
    
    Qname1 = {'EA';'ER';'ES';'ED';'EDs';'EJ';'EN';'EF';'EG';'EU';'EGstar';'EPSI';'EGmED'};
    Qname2 = 'EG';
    
    Qunit = 'h^{-1}';
    Qsr = 50*10^1;
    Qsrbar = 0.1;
    otherColor = [0.3 0.5 0.5];
     
    TEMP = load( [ PathName filesep Oname '_' Animal '_' num2str(7) 'h_olap_' num2str(0.5) filesep Oname '_backup']);
    eval(['B = TEMP.' Oname '_backup.' Qname2 ';']);
    
end


%% ContributionPlotter
if AOZ
   
    tagName = 'G';
    
    Animal = 'w140725';
    
    Oname = '';
    
    tagTensor = {'d';'do';'i'};
    globalTime = 17.25; globalOverlap = 0;
    
    Qname = {'ERcEG'; 'EScEG'; 'EDcEG' ; 'EGcEG' ; 'EAcEG' };
    
    GUI = 0;
    saveCropRegion = 0;
    
    % [X1 X2 ; Y1 Y2]
    % Qzones = {[1 1; 1 1]};
    Qzones = {[8 10; 3 5];[25 28; 2 4];[8 12; 11 13];[14 22; 9 12];[1 39; 1 19]};
    
    % Qimage = imread('\\zserver\u934\equipe_bellaiche\b_guirao\BIG_MOVIES\AOA_Outputs\WT_Contribution_MEAN_Rred_Scyan_Dgreen_sr=5-1.png');
    % Qimage = imread('C:\Users\Stephane\Desktop\FIGURES\div.png');
    Qimage = imread('X:\BIG_MOVIES\AOA_Outputs\w140725_17.25h_olap_0\Frames_circle\COA_w140725_EDscEG_d_1725_0_15h05to32h15_sr=500.png');
    
    % Qrange = {[-0.15 0.20];[-0.15 0.20];[-0.15 0.20]};
    Qrange = {};
    
    otherColor = [0.3 0.5 0.5];
    Qunit = 'h^{-1}';
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADVANCE STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% modify at your hown risk %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Format, resolution, and extention
print_format = ['-d' PLOT.extension];
print_resolution = [ '-r' num2str(PLOT.resolution)];
image_extension_output = ['.' PLOT.extension];

%% AOA - Advance
if AOA
    disp('Load Animals AIA Info for each animal ...')
            
    %%% Values to be meaned
    % NB: possibility to inverser problematique, values to not mean
    %     advantage to be more generic to any future Process. beware of specific cases such as mCD
    % TODO: Update the list with the new value name from Boris
    AOS_list = { 'AreaRatios' ; 'Rho' ; 'I'      ; 'M'       ; 'V'      };
    SM_list  = { 'AreaRatios' ; 'S'   ; 'SP'     ; 'ST'      ; 'P'      };
    VM_list  = { 'AreaRatios' ; 'U'   ; 'gradUS' ; 'gradUAS' ; 'Ianais' };
    TA_list  = { 'AreaRatios' ; 'EG'  ; 'ES'     ; 'ER'      ; 'EDs'    ; 'ED' ;  'EA' ; 'EN' ; 'EF' ; 'EJ' ; 'EJb' ; 'EDM' ; 'E' ; 'EPSI' ; 'Phi' ; 'errorPs' ; 'errorDnPs' ; 'RConds' };
    
    backupPathList = {}; % initialisation for list of path to the mean backup 
    
    %%% loading parameters and data for each animal using their respective AIA_info_
    mean_GRID_specs = cell( size(All_Animals,1), 1 );
    for n = 1:length(All_Animals)
        
        %%% load AIA info and AIA parameters for the animal       
        run([ AIAFolderName filesep 'AIA_info_' All_Animals{n} ]);
        
        %%% Defines Grid specs of all animal in order to compare them
        GRID = GridMaker(imageSize, boxSize, xyStart, gridSize, gridColor, gridLineWidth, gridOverlap); % added "gridOverlap"
        xywh = GRID.xywh;
        nx   = GRID.size(2);
        ny   = GRID.size(1);
        gridOverlap = GRID.overlap; % has been set to 0 by "GridMaker" if was empty in AIA_info
        
        olap_tag = ''; % default value: not displaying overlap in folder name when empty or 0 (for compatibility)
        if gridOverlap > 0
            olap_tag =  ['_olap_' num2str(gridOverlap)];
        end
        
        % Defines grid specific subfolder name
        gridSpecsName = ['Grid_xy_' num2str(xywh(1)) '_' num2str(xywh(2)) '_wh_' num2str(xywh(3)) 'x' num2str(xywh(4)) '_nynx_' num2str(ny) 'x' num2str(nx) olap_tag];
        
        % save the grid in the mean grid structure
        mean_GRID_specs{n} = GRID;
        
        % If it is an animal to be processed in the mean
        if ismember( All_Animals{n}, Animals )
            
            % This is an animal that we want to process, we need to load the corresponding Process Backup
            % They should be located at:
            % Animal_Folder / AIA_Info / Process_Folder / Grid_Folder / Average_Folder / Backup_Name
            %
            % This path should be stored in AIA_Parameters but we are not
            % sure that the specified average time or other information are
            % the same, so we need to rebuild the full path to find the
            % backup
            
            
            %%% Defines the "Process_Folder", "Grid_folder", "Average folder"
            %--------------------------------------------------------------
            % Process
            [~,cpt] = ismember( All_Animals{n}, Animals ); % for the AOS multichannel 
            backupPath = eval(['pathFolder_' Pname]);
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % Grid
            tag = '';
            if strcmp(Pname,'SM')
                tag = ['_' code2run '_dmu=' num2str(muAccuracy)]; % specific case for SM (for compatibility)
            end
            backupPath = [backupPath filesep gridSpecsName];
            gridLBackupFolder = [backupPath filesep 'Backups' tag];  % path to the grid backups containing the Lgrid information
            Grid_backup_rootfilename = eval(['filename_' Pname]);      % 1.9
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % Average
            % if AOA_time start and stop are define use them, else use time defined in AOT of AIA_parameters
            averageFolderName = ['Average_' num2str(animalTimeWidth) 'h_' multiTimeStart '_to_' multiTimeStop '_olap_' num2str(animalTimeOverlap)]; % moved down here (1.9)
            backupPath = [backupPath filesep averageFolderName];
            %--------------------------------------------------------------
            
            %%% Loads "mean" backup
            % NB: multiple name template possible depending on the backup:
            % mean_PROCESS_ANIMAL.mat | mean_PROCESS_ANIMAL_Backup.mat | mean_GRIDTYPE_PROCESS_ANIMAL.mat
            fullpath_backup     = [ backupPath filesep 'mean_' Pname '_' All_Animals{n} '_Backup.mat'];
            fullpath_lagrangian = [ backupPath filesep 'mean_' Pname '_LGrid_' All_Animals{n} '.mat'];
            fullpath_eulerian   = [ backupPath filesep 'mean_' Pname '_EGrid_' All_Animals{n} '.mat'];
            fullpath            = [ backupPath filesep 'mean_' Pname '_' All_Animals{n} '.mat'];
            if exist(fullpath_backup,'file')
                load(fullpath_backup)
                backupPathList = [backupPathList; fullpath_backup];
            elseif exist(fullpath,'file')
                load(fullpath)
                backupPathList = [backupPathList; fullpath];
            elseif exist(fullpath_lagrangian,'file')
                load(fullpath_lagrangian)
                backupPathList = [backupPathList; fullpath_lagrangian];
            elseif exist(fullpath_eulerian,'file')
                load(fullpath_eulerian)
                backupPathList = [backupPathList; fullpath_eulerian];
            else
                disp(['Error: backup "' fullpath '" was not found and was skipped.'])
                return;
            end
            
            %%% Store AOS multichannel name, to be used later (compatibility)
            if strcmp(Pname,'AOS')
                if length( filenameRaw ) > 1
                    idx = find( ismember( filenameRaw{1}, '_' ) );
                    for r = 1:length( filenameRaw )  % for each signal 1 a n
                        AOS_tmp_name{cpt,r} = filenameRaw{r}(1:idx-1); % get the signal name eg {'cad1_' ; 'sqh1_'} => cad1 or sqh1
                    end
                    AOS_dif_name{cpt} = intersect( AOS_tmp_name{cpt,1}, AOS_tmp_name{cpt,end}, 'stable' );
                else
                    AOS_sub_name = filenameRaw{1};
                end
            end
            
        end  % end if animal is to be process
    end  % end for each Animals
    
    
    %%% Complete the AOS list in the case of multiple signal
    if strcmp(Pname,'AOS')
        if exist('AOS_tmp_name','var')
            
            n_raw = length( filenameRaw );
            Filename_Raw_mod = cell(n_raw,1);
            ind_raw_polarity = find(polarityRawImages); % index of signal, [1 2 .... n]
            if length(Animals) == 1
                AOS_sub_name_temp = { intersect( AOS_tmp_name{1,1}, AOS_tmp_name{1,2}, 'stable' ) ; ...
                    intersect( AOS_tmp_name{1,1}, AOS_tmp_name{1,2}, 'stable' ) };
                [~,index1] = ismember( AOS_sub_name_temp{1}, AOS_tmp_name{1,1} );
                [~,index2] = ismember( AOS_sub_name_temp{2}, AOS_tmp_name{1,2} );
                AOS_sub_name = { AOS_tmp_name{1,1}(1:index1-1) ; AOS_tmp_name{1,2}(1:index2-1) };
                AOS_list = [AOS_list ; ...  % completes AOS_list with "CD_"signal names
                    ['CD_' AOS_sub_name{1} ]; ['CD_' AOS_sub_name{2} ]];
            else
                AOS_sub_name = { intersect( AOS_tmp_name{1,1}, AOS_tmp_name{2,1}, 'stable' ) ; ...
                    intersect( AOS_tmp_name{1,2}, AOS_tmp_name{2,2}, 'stable' ) };
                AOS_list = [ AOS_list ; ...  % completes AOS_list with "CD_"signal names
                    ['CD_' AOS_sub_name{1} ]; ['CD_' AOS_sub_name{2} ]];
            end
            
        end
    end
    
    %%% Loads macrocaetes positions
    if ~isempty(macrocaete_folder)
        load(macrocaete_folder);
        nbMacro = size(Landmarks_output_R,1)-1;
        macroCoord = NaN(size(Animals,1),2,nbMacro);
        macrocaetes = NaN(5,nbMacro);
        for a=1:size(Animals,1)
            idx = strmatch(Animals{a}, char(Landmarks_output_R{1,:} ));
            for m=1:nbMacro
                macroCoord(a,1,m) = Landmarks_output_R{m+1,idx}(1);
                macroCoord(a,2,m) = Landmarks_output_R{m+1,idx}(2);
            end
        end
        for m=1:nbMacro
            avg = nanmean(macroCoord(:,:,m),1);
            var = nanstd(macroCoord(:,:,m));
            if length(var) ~= 1
                r = atan2 (max( var(1), var(2)) , min( var(1), var(2)) );
                r = rad2deg(r);
                macrocaetes(:,m) = [avg(1) avg(2) var(1) var(2) r];
            else
                r = 1;
                macrocaetes(:,m) = [avg(1) avg(2) var(1) var(1) r];
            end
            
        end
        DISPLAY.macrocaetes = macrocaetes;
    end
    
    %%% Define specific possible quantities, colors, units, and path
    % Output Path
    OutputPathName = [PathName filesep 'AOA_' genericName '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) ];
    
    
            
    %% Run main process
    disp('Start AOA Processing ...')
    
    Average_Over_Animals
    
    disp('Done ... !')
end


%% DOA - Advance
if DOA
    
    disp('Load corresponding backup ...')
    
    % output folder path
    OutputFolderName = [PathName filesep 'DOA_' DOAname '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) ];
    
    LoadedBackup = cell(1,2);
    % animal 1
    disp([PathName filesep Oname '_' Animal{1} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep Oname '_backup']);
    backup = load( [PathName filesep Oname '_' Animal{1} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep Oname '_backup'] );
    LoadedBackup{1} = eval(['backup.' Oname '_backup']);
    % animal 2
    disp([PathName filesep Oname '_' Animal{2} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep Oname '_backup']);
    backup = load( [PathName filesep Oname '_' Animal{2} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep Oname '_backup'] );
    LoadedBackup{2} = eval(['backup.' Oname '_backup']);
    
    if PLOT.error
        disp(['We are calculating an error, we load the std from the first input: ' Animal{1}]);
        load([PathName filesep Oname '_' Animal{1} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep Oname '_stdMAP']);
    end
    
    disp('Start DOA processing ...')
    
    Difference_Over_Animals;
    
    disp('Done ... !')
end


%% COQ - Advance
if COQ
    
    
    
    image_extension_output = ['.' PLOT.extension];
    print_format = ['-d' PLOT.extension];
    print_resolution = ['-r' num2str(PLOT.resolution)];
    
    
    
    InputPathName  = [PathName filesep Oname '_' Animal '_' num2str(COQ_time_width) 'h_olap_' num2str(COQ_time_overlap) ];
    OutputPathName = [InputPathName ];
    
    disp('Loading corresponding backup ...')
    load( [InputPathName filesep Oname '_backup'] );
    eval( ['BACKUP = ' Oname '_backup;']);
    
    disp('Start Processing ...')
    Calcul_Over_Quantities;
    
    disp('Done ... !')
end


%% COA - Advance
if COA
    
    InputPathName  = [ PathName filesep Oname '_' Animal '_' num2str(COA_time_width) 'h_olap_' num2str(COA_time_overlap) ];
    OutputPathName = [InputPathName ];
    
    if deviator && isotropic
        deviator = 0;
        isotropic = 0;
    end
    
    tensorTag = '';
    if deviator
        tensorTag = '_d';
    elseif isotropic
        tensorTag = '_i';
    end
    denominatorTag = '';
    if strcmp(denominator,'norm2')
        denominatorTag = '_n2';
    elseif strcmp(denominator,'norm2W')
        denominatorTag = '_n2W';
    elseif strcmp(denominator,'mean')
        denominatorTag = '_m';
    end
    
    disp('Loading corresponding backup ...')
    load( [InputPathName filesep Oname '_backup'] );
    eval(['BACKUP = ' Oname '_backup;']);
    
    disp('Start COA processing ...')
    
    Contraction_Over_Animal_new;
    
    disp('Done ... !')
end

%% AOZ - Advance