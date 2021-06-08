% AIA_MultiOperation
%
% Script containing all necessary info required to run Operation over Animal
% and is aimed to extract most of those information respectively from the 
% different AIA_Info of each animals.
%
% definitions and conventions
% Quantities : EG, ED, EA, I, V, S, ... ; shortened into Q
% Processes  : AOS, TA, SM, ... ;         shortened into P
% Operations : AOA, DOA, ... ;            shortened into O
%
% With exception of AOA and DOA, all operations run on single and average
% movies. To run a single movie, one must leave the Oname empty.
% Single animal operation can only be done on quantities from the same
% process (e.g. cannot do EG and S on the same run)
%
% ACHTUNG! POA can involve 2 different animal. The corresponding backup and
% plot will be named by the secondary animal, but saved in the primary
% animal folder.
%
% ACHTUNG again! COQ on single animal will directly update the raw backup
% generated by Boris framework
%
% TODO:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% * take advantage of the fact that now all quantities from all processes TA,SM,AOS... are gathered in the SAME backup!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% * finishing integrating COQ quantities in the other operations
% * GUI interface for selecting region in LTA
% * normalisation of tensor norm in LTA (needed?)
% * Clean and add Tensor Correlator
% * Test all the possible options (should work fine ... I guess)
%
% Version 2.1
% Stephane Rigaud
% Boris Guirao


clear all; close all; clc;
AIA_call = 1;               % Flag preventing AIA_parameters to run other programs (DO NOT MODIFIY)


%% Operation Name (Oname) Selection %%

% Boolean in order to select which Operation will be run ;
% NB: /!\ only run one at a time, otherwise the planete will implode

AOA = 0;               % Average Over Animal (single animal only)
DOA = 0;               % Difference Over Animal (averaged animal only)
COQ = 0;               % Calcul Over Quantities
POA = 0;               % Projection Over Animal
LTA = 1;               % Local Time Analysis

TC  = 0;               % Tensor Corelator (in progress)

% Global output path: all outputs will be saved there in a structure file.
PathName = 'D:\BigMovies\AOA_WT_stress_fibers';
AIAFolderName = 'C:\Users\Boris\ownCloud\Matlab\AIA_infos'; % path to the folder containing the AIA_info_X


%% Main Common Parameters %%

Animal = 'meanTRBL_grid26h';        % name of average OR actual animal name (BIG1, TRBL4...) when processing single animal
% Animal ='meanWT_grid26h';
% Animal = 'BIGwt2';
%     Animal = 'R23';
% Animal = 'TRBL8';

singleAnimal = false;               % processing of single animal or average ?

% Quantities to be PLOTTED (ALL "Pname" quantities are CALCULATED) in AOA, DOA and PROJECTED (calculated and plotted) in POA
Qname = {'EA';'ER';'ES';'ED';'EG'};  % Quantity names
% Qname = {'EA';'ER';'ES';'ED';'EJ';'EN';'EF';'EG'};  % Quantity names
%     Qname = {'M' 'Rho' };  


%% Time information (mod 2.1) %%

% animalTimeWidth   = 20;
% animalTimeOverlap = 0;  % to make a plot at every single frame
animalTimeWidth   = 2;
animalTimeOverlap = 0.96;  % to make a plot at every single frame

% multiTimeStart    = '00h00';
% multiTimeStop     = '20h00';
multiTimeStart    = '12h00';
multiTimeStop     = '32h00'; % 28 for wt3

gridTime = '26h00';     % time APF used to draw grid and determine cell patches to track (2.1) 
delta_t = 5;            % LTA (moved here 2.1)
temperature = 29;       % to correct "delta_t" (2.1)
timeShift = 2;      % in HOURS, correction of reference time (2.1)


%% Projection (moved 2.1) %%

% NB: relevant for POA, LTA

% if projection is involved, specify ONTO what here below (1st run POA, then LTA):
uAnimal = 'meanTRBL_grid26h';
% uAnimal = 'TRBL8';
% uAnimal = 'R23';
% uAnimal ='meanWT_grid26h';
% uAnimal = 'BIGwt2';

uTimeWidth = 20;
uTimeOverlap = 0;
uQname = 'EG';      % quantity used for projection (specific to POA)
uOname = 'AOA';      % Operation name (AOA,DOA...) that led to "uAnimal" (overridden by "empty" when single animal processing"
% NB: if "uAnimal" and "uOname" are empty, we use the OnamePOA and Animal from the projection value


%% Common Plot parameters for tensor maps (not LTA) %%

PLOT.plot = 0;                      % plot the selected quantities (Qname)
PLOT.print = 1;                     % print the plot
PLOT.significance = 0;              % plot significance
PLOT.SignOpacityMap = [0.5 0.5];    % significance opacity [iso dev]
PLOT.extension = 'png';             % plot image type (png | svg)
PLOT.resolution = 300;              % plot resolution
PLOT.boxSize = [122 128/2]/0.322;   % JESUS
% PLOT.boxSize = [256 256];        % box size of the average grid
PLOT.macrocaetes = 1;               % plot macrocaetes position

QplotType = 'split+';               % split+/-, dev+/-, circle, merge (moved 2.1)


%% AOA - Average over Animals %%

% This process calculate the average of N individual (N can be 1) together.
% The spatial and temporal average should be the same between the animal
% in order to average them together (see AOT)

% Specify the list of all animals name that will be used in the generation of the commong average grid
% NB: The list should contain ALL animals that will be compared and not only the animal to be ploted
% by default: provide all the possible animals


% ALL animals to consider to build commong grid:
%------------------------------------------------------------
% allAnimals = { 'TRBL1' ; 'TRBL4' ; 'TRBL7' ; 'TRBL9' ; 'TRBL8l'  ; 'TRBL8r'  }; %#ok<*UNRCH> % animal name list
%     allAnimals = {'wt3'};
allAnimals = { 'BIG_1' ; 'BIG_5' ; 'BIG_6' ; 'BIG_N5m' ; 'BIGwt2r' ; 'BIGwt2l'};
%------------------------------------------------------------


% Animals to take into account in the average:
%------------------------------------------------------------
% avgAnimals = { 'TRBL1' ;'TRBL4' ; 'TRBL7' ; 'TRBL9' ; 'TRBL8l' ; 'TRBL8r'};
%    avgAnimals = {'wt3'};
avgAnimals = { 'BIG_1' ; 'BIG_5' ; 'BIG_6' ; 'BIG_N5m' ; 'BIGwt2r' ; 'BIGwt2l'};


% ESCARGOT FULL 12 APF
% allAnimals = {'esgGfp_12apf_1';'esgGfp_12apf_2';'esgGfp_12apf_3';'esgGfp_12apf_4'};
% avgAnimals = {'esgGfp_12apf_1';'esgGfp_12apf_2';'esgGfp_12apf_3';'esgGfp_12apf_4'};


% ESCARGOT FULL 24 APF
% allAnimals = {'esgGfp_24apf_1';'esgGfp_24apf_2';'esgGfp_24apf_3'};
% avgAnimals = {'esgGfp_24apf_1';'esgGfp_24apf_2';'esgGfp_24apf_3'};
%------------------------------------------------------------

PnameAOA = 'TA'; % Process name: SM,VM,AOS,TA ** TO BE REMOVED **
%     PnameAOA = 'AOS';         % Process name: SM,VM,AOS,TA



%% DOA - Difference over Animals

% This process calculate the difference between two individual.
% In this case, an individual is the output of a Process (AOA, DOA, COQ, other)
% In the current state, both individual MUST BE COMING from the same Operation (AOA - AOA, DOA - DOA, but no AOA - DOA)
% As this is a simple A - B operation, both individual must have the same spatial and temporal size, with same overlap

% NB: to perform the difference between a quantity Q from 2 single animals (Qa1-Qa2), one HAS TO first "average" each of
% them through AOA in order to build a grid COMMON to BOTH animals.


OnameDOA = 'AOA'; % Operation name (single animal don't work)

% Animal 1 minus Animal 2: the name of the output animal will be "A1-A2"
deltaAnimals  = { 'meanWT3' ; 'meanWT3int'};

PnameDOA = 'TA';            % Process name ** TO BE REMOVED **

PLOT.Qparameters = 1;
% if Qparameters is off, the following parameters are used for the plot
Qcolor = [0.3 0.5 0.5];
Qunits = 'h^{-1}';
Qsr = 6e3;
Qsb = 1e-2;
%     Qsr = 4e3;
%     Qsrbar = 2e-2;
QKillTr = 0;


%% COQ - Calcul over Quantities

% COQ allow to perform simple mathematical operation, such as substraction or addition, on quantities or an animal.
% The difference with DOA is that all operation are performed between quantities from the SAME animal, e.g. A.EG - A.ED
% The results is saved in the animal backup as a new quantity
  
OnameCOQ = 'AOA'; % Operation name (overridden by "empty" when single animal processing)

% Quantities involved in operation A - B, A + B, ...
% NB: A and B should be from the same process (TA, AOS...) AND have same units
QnameA = 'EG';
QnameB = 'ED';

CalculType = 'minus'; % choose between 'minus', 'plus'


%% POA - Projection over Animals %%

% POA projects the quantity A on the unitary tensor of the quantity B
% A and B can be from different animal with different temporal averaging
% However their grid must be identical

% type of projection to be plotted in maps (only one at a time)
tagPlot = 'd';  % choose between 'd', 'i', 'do' for, respectively deviator, isotropic, orthogonal deviator

% project this:
OnamePOA = 'AOA';              % Operation name (AOA,DOA...) (overridden by "empty" when single animal processing)

QunitsPOA = 'h^{-1}';
QsrPOA    = 50*10^1;
QsbPOA    = 0.1;



%% LTA - Local Time Analysis %%


OnameLTA = 'POA';         % Operation name (AOA, DOA, POA...)

QnameLTA = {'EGdotEG','ERdotEG','ESdotEG','EDdotEG','EAdotEG'}; % list of quantities to be PLOTTED
%     Qname = {'EG','ED','ER','ES','EA','M','I','Rho'};
%     Qname = {'EDdotEG','ERdotEG','ESdotEG','EAdotEG'};
%     Qname = {'ID2'};

% range of zone xy coordinates (not ij!!!) to be included (those coordinates start at [1 1]):
%     Qareas = {[2 2; 2 2]}; % {[xMin xMax ; yMin yMax]} of box coodinates % TRBL8
Qareas = {[2 2; 1 1]}; % {[xMin xMax ; yMin yMax]} of box coodinates % meanTRBL & meanWT
%     Qareas = {[1 1; 1 1]}; % {[xMin xMax ; yMin yMax]} of box coodinates % R23
%     Qareas = {[2 2; 3 3]}; % {[xMin xMax ; yMin yMax]} of box coodinates % BIGwt2
% NB: Examples: {[1 6; 1 5]} will process all boxes from (1,1) to box (6,5) (in xy coordinates)
% NB: Examples: {[6 6; 5 5]} will ONLY process the box (6,5) (in xy coordinates)

%     CropRegionGrid = Qareas{z};
%     Quantity(CropRegionGrid(2,1):CropRegionGrid(2,2),CropRegionGrid(1,1):CropRegionGrid(1,2),:,:);
%     filename = [ '(' num2str(CropRegionGrid(1,1)) ',' num2str(CropRegionGrid(2,1)) ')(' num2str(CropRegionGrid(1,2)) ',' num2str(CropRegionGrid(2,2)) ')'];

% Axis range, leave empty for automatic range selection, otherwise like this: {[400 410]}
%     Qrange = {};
Qrange = {0.07*[-1 1]}; % BIG & TRBL
%     Qrange = {0.12*[-1 1]};     % R23

Normalisation = ''; % normalisation over region, animal, or no normalisation (none)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADVANCED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% modify at your hown risk %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Overriding values by "empty" when processing single animal (2.1)
if singleAnimal
    OnameCOQ = '';
    OnamePOA = '';
    uOname = '';
end


%%% Time Rescaling: correcting "delta_t" according to temperature (2.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if temperature == 29
    delta_t = delta_t / 0.9;          % to compare velocities, constriction rates... between 25 and 29 movies, dt must be corrected
    disp(['WARNING: "AIA_info" parameter "delta_t" has been updated to ' num2str(delta_t) ' since temperature = 29!'])
elseif temperature ~= 25
    disp(['ERROR: "AIA_info" parameter "temperature" (here ' num2str(temperature) ') can only be 25 or 29!'])
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Format, resolution, and extention
printFormat     = ['-d' PLOT.extension];
printResolution = ['-r' num2str(PLOT.resolution)];
imageExtension  = ['.'  PLOT.extension];

% Store time information in struct more simple managment
TIME.animalTimeOverlap = animalTimeOverlap;
TIME.animalTimeWidth   = animalTimeWidth;
TIME.multiTimeStart    = multiTimeStart;
TIME.multiTimeStop     = multiTimeStop;

CustomColors;       % defines usual set of colors
AllQsColorsUnits;   % associate quantities with specific colors

% display parameters to be used during plot
minimalInfoDisplay = false;
% minAEV = 0.0001;
fontSize = 20;
EVstyles = {'-' '--'};       % ONLY relevant for "merged" display type: styles to display ellipse axes representing tensor eigenvalues (default {'-' ':'})
pointSize = 2;
signOpacities = [0.7 0.3];   % ONLY relevant for "split+/-" and "circle" display types: specifies opacity of positive(white) and negative(black) disks, respectively.
lineWidth = 1.5;             % for circle, bars and ellipses (1.5 ok with BIG movies)
gridDisplay = false;         % Lagrangian grid ALWAYS displayed (6.0)
gridColor = black;
gridLineWidth = 0.5;        % only matters for Egrid, Lgrid thickness specified in LGridPlotter
imageFading = 0.6;
scaleBarWidth = 1;

% Scale values according to the information to be ploted
sr_AOS        = { 4e2 ; 50 ;0.8*[1 3]; 100 ; [4 30] };   % sets ratio setting size of ellipses or bars in tensor representation for M,I. % 2 hours
% sr_AOS        = { 5e2 ; 50 ;[1 3]; 100 ; [4 30] };   % sets ratio setting size of ellipses or bars in tensor representation for M,I. % 14h/20h hours
srbar_AOS     = [ 0.1 ;  2 ; 50 ;  1  ;    10  ];   % scale bar lengths for each contribution
killtrace_AOS = [  0  ;  0 ;  0 ;  1  ;     1  ];   % Will set average compartment trace to 0 in the plots (mean isotropic part = 0). Choose this when tensors are known up to an additive constant
% contributions   Rho    I    M    V       CD

sr_SM        = { 1500 ; 1000 ; 1000 ; 1000 };
srbar_SM     = [  0.1 ; 0.05 ; 0.05 ; 0.05 ];
killtrace_SM = [    1 ;    1 ;    1 ;    1 ];
% contributions     S     SP     ST      P

% eLife paper:
%-----------------------------------------------------------------------------------------------------------------------------------------
% sr_TA        = {  4e3 ;  2e3 ;  4e3 ;  8e3 ;  2e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 };  % ISO average 14h/20h
% sr_TA        = {  4e3 ;  4e3 ;  4e3 ;  8e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 };  % DEV average 14h/20h
% srbar_TA     = [ 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ; 2e-2 ];
%-----------------------------------------------------------------------------------------------------------------------------------------

%%% 14h/20h averages:
%-----------------------------------------------------------------------------------------------------------------------------------------
% sr_TA        = {  4e3 ;  1.5e3 ;  4e3 ;  8e3 ;  1.5e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 };  % ISO average 14h/20h
sr_TA        = 6e3; sr_TA = num2cell(sr_TA*ones(15,1));  % DEV average 14h
% sr_TA        = 4e3; sr_TA = num2cell(sr_TA*ones(15,1));  % DEV average 14h
srbar_TA     = 2e-2*[ 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ];
%-----------------------------------------------------------------------------------------------------------------------------------------
%%% 2h averages:
%-----------------------------------------------------------------------------------------------------------------------------------------
% sr_TA        = {  6e3 ;  [1.5e3 6e3] ;  6e3 ;  6e3 ;  [1.5e3 6e3] ;  6e3 ;  6e3 ;  6e3 ;  6e3 ;  6e3 ;  6e3 ;  6e3 ;  6e3 ;  6e3 ;  6e3 };  % ISO&DEV average 2h
% srbar_TA     = 2e-2*[ 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ];
% STEPHANE:
% sr_TA        = 2e3; sr_TA = num2cell(sr_TA*ones(15,1));  % DEV average 2h
% sr_TA        = {  2e3 ;  0.5e3 ;  2e3 ;  2e3 ;  0.5e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 };  % ISO average 2h
% sr_TA        = {  2e3 ;  [0.5e3 2e3] ;  2e3 ;  2e3 ;  [0.5e3 2e3] ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 ;  2e3 };  % ISO&DEV average 2h
% srbar_TA     = 5e-2*[ 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ; 1 ];
%-----------------------------------------------------------------------------------------------------------------------------------------

killtrace_TA = [    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ;    0 ];
% contributions     G      S      R     Ds      D      A      N      F      J     Jb     DM      U     G*    PSI   PhiU

sr_VM        = { 10 ; 2000 ; 2000 };
srbar_VM     = [  1 ; 0.02 ; 0.02 ];
killtrace_VM = [  0 ;    0 ;    0 ];
% contributions   U  gradUS  gradUAS

sr_GEP        = {  150 };
srbar_GEP     = [  0.5 ];
killtrace_GEP = [  0 ];
% contributions   ID

%% AOA - Advanced %%

if AOA
    disp('Load Animals AIA Info for each animal ...')
    
    % Values to be meaned
    % NB: possibility to inverser problematique, values to not mean
    %     advantage to be more generic to any future Process. beware of specific cases such as mCD
    % TODO: Update the list with the new value name from Boris
    %AOS_list = { 'AreaRatios' ; 'Rho' ; 'I'      ; 'M'       ; 'V'      };
    %SM_list  = { 'AreaRatios' ; 'S'   ; 'SP'     ; 'ST'      ; 'P'      };
    %VM_list  = { 'AreaRatios' ; 'U'   ; 'gradUS' ; 'gradUAS' ; 'Ianais' };
    %TA_list  = { 'AreaRatios' ; 'EG'  ; 'ES'     ; 'ER'      ; 'EDs'    ; 'ED' ;  'EA' ; 'EN' ; 'EF' ; 'EJ' ; 'EJb' ; 'EDM' ; 'E' ; 'EPSI' ; 'Phi' ; 'errorPs' ; 'errorDnPs' ; 'RConds' };
    rejectList = {'xywh';'size';'overlap';'centroids';'ULCs'; 'fullImage';'color';'lineWidth';'coordinates'; 'TimeArray';'FrameArray'};
    
    backupPathList = {}; % initialisation for list of path to the mean backup
    macroCoord = NaN(14,2,length(allAnimals));
    
    % loading parameters and data for each animal using their respective AIA_info_
    mean_GRID_specs = cell( size(allAnimals,1), 1 );
    for n = 1:length(allAnimals)
        
        disp(allAnimals{n});
        
        % If it is an animal to be processed in the mean
        if ismember( allAnimals{n}, avgAnimals )
            [~,cpt] = ismember( allAnimals{n}, avgAnimals ); % for the AOS multichannel
            [GRID backupPath AOSvariable, RAW] = SingleAnimalLoader(AIAFolderName, allAnimals{n}, PnameAOA, TIME, true );
            
            % Store backup path to be loaded later
            backupPathList = [backupPathList backupPath];
            
            % THIS IS NOT WORKING -----------------------------------------
            % Store AOS multichannel name, to be used later (compatibility)
%             if strcmp(PnameAOA,'AOS') 
%                 AOS_tmp_name{cpt,:} = AOSvariable.AOS_tmp_name;
%                 AOS_dif_name{cpt} = AOSvariable.AOS_dif_name;
%                 AOS_sub_name = AOSvariable.AOS_sub_name;
%             end
            % THIS IS NOT WORKING -----------------------------------------
            
%             eval(['GRID_' PnameAOA '_' allAnimals{n} '= load(backupPath);']);
            load(backupPath)
            
            %%% load macrocaete position
            %--------------------------------------------------------------
            nbMacro = 8;
            if strcmp(RAW.halfNotum,'b')
                nbMacro = 16; % upped from 14 (2.1)
            end
            macroCoord(1:nbMacro,:,n) = LandmarkLoader(RAW.rawPathFolder, allAnimals{n}, RAW.halfNotum, gridTime); % added "gridTime" (2.1)
%             macroCoord(1:nbMacro,:,n) = LandmarkLoader(RAW.rawPathFolder, allAnimals{n}, RAW.halfNotum);
            coefXY = PLOT.boxSize ./ RAW.boxSize;
            macroCoord(:,1,n) = (macroCoord(:,1,n) - RAW.xyStart(1)) .* coefXY(1);
            macroCoord(:,2,n) = (macroCoord(:,2,n) - RAW.xyStart(2)) .* coefXY(2);
            %--------------------------------------------------------------
        else
            [GRID, ~, ~, ~] = SingleAnimalLoader(AIAFolderName, allAnimals{n}, PnameAOA, TIME, false );
        end
        
        % save the grid in the mean grid structure
        mean_GRID_specs{n} = GRID;
        gridOverlap = GRID.overlap;
    end
    
    % Complete the AOS list in the case of multiple signal
%     if strcmp(PnameAOA,'AOS') % THIS IS NOT WORKING
%         if exist('AOS_tmp_name','var')
%             
%             n_raw = length( filenameRaw );
%             Filename_Raw_mod = cell(n_raw,1);
%             ind_raw_polarity = find(polarityRawImages); % index of signal, [1 2 .... n]
%             if length(avgAnimals) == 1
%                 AOS_sub_name_temp = { intersect( AOS_tmp_name{1,1}, AOS_tmp_name{1,2}, 'stable' ) ; ...
%                     intersect( AOS_tmp_name{1,1}, AOS_tmp_name{1,2}, 'stable' ) };
%                 [~,index1] = ismember( AOS_sub_name_temp{1}, AOS_tmp_name{1,1} );
%                 [~,index2] = ismember( AOS_sub_name_temp{2}, AOS_tmp_name{1,2} );
%                 AOS_sub_name = { AOS_tmp_name{1,1}(1:index1-1) ; AOS_tmp_name{1,2}(1:index2-1) };
%                 AOS_list = [AOS_list ; ...  % completes AOS_list with "CD_"signal names
%                     ['CD_' AOS_sub_name{1} ]; ['CD_' AOS_sub_name{2} ]];
%             else
%                 AOS_sub_name = { intersect( AOS_tmp_name{1,1}, AOS_tmp_name{2,1}, 'stable' ) ; ...
%                     intersect( AOS_tmp_name{1,2}, AOS_tmp_name{2,2}, 'stable' ) };
%                 AOS_list = [ AOS_list ; ...  % completes AOS_list with "CD_"signal names
%                     ['CD_' AOS_sub_name{1} ]; ['CD_' AOS_sub_name{2} ]];
%             end
%             
%         end
%     end
    
    % Loads macrocaetes positions
%    if ~isempty(macrocaete_folder)
%         load(macrocaete_folder);
%         nbMacro = size(Landmarks_output_R,1)-1;
%         macroCoord = NaN(size(avgAnimals,1),2,nbMacro);
%         macrocaetes = NaN(5,nbMacro);
%         for a=1:size(avgAnimals,1)
%             idx = strmatch(avgAnimals{a}, char(Landmarks_output_R{1,:} ));
%             for m=1:nbMacro
%                 macroCoord(a,1,m) = Landmarks_output_R{m+1,idx}(1);
%                 macroCoord(a,2,m) = Landmarks_output_R{m+1,idx}(2);
%             end
%         end
%         for m=1:nbMacro
%             avg = nanmean(macroCoord(:,:,m),1);
%             var = nanstd(macroCoord(:,:,m));
%             if length(var) ~= 1
%                 r = atan2 (max( var(1), var(2)) , min( var(1), var(2)) );
%                 r = rad2deg(r);
%                 macrocaetes(:,m) = [avg(1) avg(2) var(1) var(2) r];
%             else
%                 r = 1;
%                 macrocaetes(:,m) = [avg(1) avg(2) var(1) var(1) r];
%             end
%             
%         end
        macrocaetes = NaN(5,14);
        macrocaetes(1:2,:) = nanmean(macroCoord,3)';
        DISPLAY.macrocaetes = macrocaetes;
    %end
    
    % Output Path
    OutputPathName = [PathName filesep 'AOA_' Animal '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) ];
    
    % Run main process
    disp('Start AOA Processing ...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AverageOverAnimals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Done ... !')
end


%% DOA - Advanced %%

if DOA
    
    DOAname = [deltaAnimals{1} '-' deltaAnimals{2}];
    
    disp('Load corresponding backup ...')
    
    LoadedBackup = cell(1,2);
    if isempty(OnameDOA)
        
        disp('ERROR: DOA cannot work on raw animal because of grid correspondances.')
        disp('TIPS : Do a fake average of the animal for fixing the grid issue.')
        
        %         % animal 1
        %         [GRID backupPathList AOSvariable] = SingleAnimalLoader(AIAFolderName, deltaAnimals{1}, PnameDOA, TIME, true );
        %         disp(backupPathList)
        %         load(backupPathList);
        %         LoadedBackup{1} = eval(['GRID_' PnameDOA '_' deltaAnimals{1}]);
        %
        %         % animal 2
        %         [GRID backupPathList AOSvariable] = SingleAnimalLoader(AIAFolderName, deltaAnimals{2}, PnameDOA, TIME, true );
        %         disp(backupPathList)
        %         load(backupPathList);
        %         LoadedBackup{2} = eval(['GRID_' PnameDOA '_' deltaAnimals{2}]);
        
    else
        
        % animal 1
        disp([PathName filesep OnameDOA '_' deltaAnimals{1} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep OnameDOA '_backup']);
        LoadedBackup{1} = load( [PathName filesep OnameDOA '_' deltaAnimals{1} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep OnameDOA '_backup'] );
        
        % animal 2
        disp([PathName filesep OnameDOA '_' deltaAnimals{2} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep OnameDOA '_backup']);
        LoadedBackup{2} = load( [PathName filesep OnameDOA '_' deltaAnimals{2} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep OnameDOA '_backup'] );
        
    end
    
    
    % output folder path
    OutputFolderName = [PathName filesep 'DOA_' DOAname '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) ];
    
    %     if PLOT.error %%% ????????
    %         disp(['We are calculating an error, we load the std from the first input: ' deltaAnimals{1}]);
    %         load([PathName filesep OnameDOA '_' deltaAnimals{1} '_' num2str(animalTimeWidth) 'h_olap_' num2str(animalTimeOverlap) filesep OnameDOA '_stdMAP']);
    %     end
    
    disp('Start DOA processing ...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DifferenceOverAnimals;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Done ... !')
end


%% COQ - Advanced %%

if COQ
    
    disp('Loading corresponding backup ...')
    if isempty(OnameCOQ)
        PLOT.significance = 0;
        PLOT.macrocaetes = 0;
        Pname = GetPname(QnameA);
        [~, backupPathList, ~, RAW] = SingleAnimalLoader(AIAFolderName, Animal,  'AOT',  TIME,  true );
        BACKUP = load( backupPathList );
        [OutputPath,backupPathName,~]  = fileparts(backupPathList);
        backupName = ['GRID_' Pname '_' Animal];
    else
        InputPathName  = [PathName filesep OnameCOQ '_' Animal '_' num2str(TIME.animalTimeWidth) 'h_olap_' num2str(TIME.animalTimeOverlap) ];
        BACKUP = load( [InputPathName filesep OnameCOQ '_backup'] );
        OutputPath = InputPathName ;
        backupName = [OnameCOQ '_backup']; 
        backupPathName  = [OnameCOQ '_backup'];
    end
    
    disp('Start Processing ...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CalculationOverQuantities;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Done ... !')
end


%% POA - Advanced %%

if POA
    
    tagTensor = {'i';'d';'do'};
    tagProjectionTime = [strrep(num2str(uTimeWidth),'.','') '_' strrep(num2str(uTimeOverlap),'.','')];
    
    QplotType = 'circle';
    QKillTr = 0;
    
    Pname  = GetPname(Qname{1});
    uPname = GetPname(uQname  );
    
    uTIME = TIME;
    uTIME.animalTimeWidth = uTimeWidth;
    uTIME.animalTimeOverlap = uTimeOverlap;
    
    if isempty(OnamePOA)
        
        % projection is a single animal: 
        % - Animal must be the same as uAnimal
        % - OnamePOA and uOname must be empty
        % - Qname can be different than uQname
        
        uOname = '';
        uAnimal = Animal;
        PLOT.significance = false;
        
        [~, backupPathList, ~, RAW] = SingleAnimalLoader(AIAFolderName, Animal,  'AOT',  TIME, true );
        [~, ubackupPathList, ~, ~]  = SingleAnimalLoader(AIAFolderName, uAnimal, 'AOT', uTIME, true );
        
        disp( 'Loading projection backup ...' );
        disp( backupPathList );
        BACKUP = load( backupPathList );
        
        disp( ['Loading projector backup ...' uQname] );
        disp( ubackupPathList );
        TEMPBACKUP = load( ubackupPathList );
        eval( ['B = TEMPBACKUP.' uQname ';'] );
        
        animalProjectionFolder = [Animal '_' num2str(TIME.animalTimeWidth) 'h_olap_' num2str(TIME.animalTimeOverlap)];
        
%         tempSubFolder = fileparts(backupPathList);
%         gridAnimalFolder = fileparts(tempSubFolder);
%         rawPathFolder = fileparts(fileparts(fileparts(tempSubFolder)));

    else
        
        % projection is an averaged animal: 
        % - check if uAnimal is specified, else uAnimal = Animal
        % - check if uOname is specified, else uOname  = OnamePOA
        % - Qname can be different than uQname
        
        if isempty(uAnimal), uAnimal = Animal; end
        if isempty(uOname),  uOname  = OnamePOA;  end
        
        animalProjectionFolder = [OnamePOA  '_' Animal  '_' num2str(TIME.animalTimeWidth)  'h_olap_' num2str(TIME.animalTimeOverlap) ];
        animalProjectorFolder  = [uOname '_' uAnimal '_' num2str(uTIME.animalTimeWidth) 'h_olap_' num2str(uTIME.animalTimeOverlap)];
        
        disp( 'Loading projection backup ...' )
        disp( [PathName filesep animalProjectionFolder filesep OnamePOA '_backup'] );
        BACKUP = load( [PathName filesep animalProjectionFolder filesep OnamePOA '_backup'] );
        
        disp( ['Loading projector value ...' uQname] )
        disp( [PathName filesep animalProjectorFolder filesep OnamePOA '_backup'] );
        TEMPBACKUP = load( [PathName filesep animalProjectorFolder filesep OnamePOA '_backup'] );
        eval( ['B = TEMPBACKUP.' uQname ';'] );
        
    end
    
    outputPath = [PathName filesep animalProjectionFolder];
    
    disp('Start POA processing ...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ProjectionOverAnimals;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Done ... !')
end


%% LTA - Advanced %%

if LTA
    
    
    tagTensor = {''};
    tagProjection = '';
    tagBackup = 'backup';
    if strcmp(OnameLTA,'POA')
        tagTensor = {'_d'; '_i'; '_do'};
        tagProjection = ['_' strrep(num2str(uTimeWidth),'.','') '_' strrep(num2str(uTimeOverlap),'.','')];
        tagBackup = uAnimal;
    end
    backupName = [OnameLTA '_' tagBackup tagProjection];
    backupPathName = [OnameLTA '_' tagBackup tagProjection];
    
    if singleAnimal
        
        if strcmp(OnameLTA,'POA')
            backupPathList = [PathName filesep Animal '_' num2str(TIME.animalTimeWidth) 'h_olap_' num2str(TIME.animalTimeOverlap)];
            backupPathList = [backupPathList filesep backupPathName ];
        else
            Pname = GetPname(QnameLTA{1});
            [~, backupPathList, ~, RAW] = SingleAnimalLoader(AIAFolderName, Animal,  'AOT',  TIME,  true );
            backupName = ['GRID_' Pname '_' Animal];
        end
        
        OutputPath = fileparts(backupPathList);
        
        disp( 'Loading backup ...' );
        disp( backupPathList );
        BACKUP = load( backupPathList );
          
    else
        
        if strcmp(OnameLTA,'DOA')
            InputPathName = [PathName filesep 'DOA_' Animal '_' num2str(TIME.animalTimeWidth) 'h_olap_' num2str(TIME.animalTimeOverlap)];
        else
            InputPathName = [PathName filesep 'AOA_' Animal '_' num2str(TIME.animalTimeWidth) 'h_olap_' num2str(TIME.animalTimeOverlap)];
        end
        
        backupPathList = [InputPathName filesep backupPathName];
        
        OutputPath = fileparts(backupPathList);
        
        disp( 'Loading backup ...' );
        disp( backupPathList );
        BACKUP = load( backupPathList );
        
    end
    
    disp('Start LTA processing ...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LocalTimeAnalysis;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Done ... !')
    
end


%% History %%

% 04/01/2017: 2.1 (Boris)
% - fixed "delta_t" that was NOT updated for animals at 29!! (impacted LTA plots)
% - changed all "Oname" defined in each Operation section into "OnameDOA", "OnameCOA", "OnamePOA", "OnameLTA" so Operations can run
% sequentially (without having to run each SEPARATELY)
% - defined common "Qname" defined at beginning (except for LTA)
% - added parameter "gridTime" in time parameter section (and as additional "LandmarkLoader" argument) now that grid can be defined at
% arbitrary timepoint.
% - increased "nbMacro" from 14 to 16
% - added "timeShift" to correct for rotation peak delay in time APF (when later defining "startAPF" and "endAPF")
% - QplotType = split+/-, dev+/-, circle, merge became a common parameter
% - in AOA section: "All_Animals" became "allAnimals"; "Animals" became "avgAnimals"; "genericName" became "Animal";"Pname" became "PnameAOA"
% - in DOA section: "Animal" became "deltaAnimals"; "Pname" became "PnameDOA"
% - in POA section: Qunits, Qsr, Qsb became QunitsPOA, QsrPOA, QsbPOA
% - "Animal" is now common to all operations and has been moved to the top
% - Qsrbar became Qsb
