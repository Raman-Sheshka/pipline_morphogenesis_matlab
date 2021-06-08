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
% NB1: POA can involve 2 different animals. The corresponding backup and
% plot will be named by the secondary animal, but saved in the primary
% animal folder.
%
% NB2: COQ on single animal will directly update the AOT backup generated by Boris framework
%
% TO DO:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% * take advantage of the fact that now all quantities from all processes TA,SM,AOS... are gathered in the SAME backup!
% * allow to project tensors along xx and yy directions
% * symmetrize average tensor maps % midline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% * finishing integrating COQ quantities in the other operations
% * GUI interface for selecting region in LTA
% * normalisation of tensor norm in LTA (needed?)
% * Clean and add Tensor Correlator
% * Test all the possible options (should work fine ... I guess)
%
% Version 2.3
% Stephane Rigaud
% Boris Guirao


clear all; close all; clc;
AIA_call = 1;               % Flag preventing AIA_parameters to run other programs (DO NOT MODIFIY)


%% Operation Name (Oname) Selection %%

% Boolean in order to select which Operation will be run ;
% NB: /!\ only run one at a time, otherwise the planete will implode

AOA = 0;               % Average Over Animal (single animal only)
DOA = 0;               % Difference Over Animal (averaged animal only)  (NOT up to date for AOT folders)
COQ = 0;               % Calcul Over Quantities                         (NOT up to date for AOT folders)
POA = 0;               % Projection Over Animal
LTA = 1;               % Local Time Analysis 

TC  = 0;               % Tensor Corelator (in progress)

% Global output path: all outputs will be saved there in a structure file.
% PathName = 'D:\BigMovies\AOA_WT_128grid_olap0.5';
PathName = 'D:\BigMovies\AOA_128grid_olap0.5';
AIAFolderName = 'C:\Users\Boris\ownCloud\Matlab\AIA_infos'; % path to the folder containing the AIA_info_X


%% Main Common Parameters %%

% Animal = 'meanTRBL_grid26h';        % name of average OR actual animal name (BIG1, TRBL4...) when processing single animal
Animal ='meanWT';
% Animal ='meanWT_grid26h';
% Animal = 'BIGwt2';
% Animal = 'R23';
% Animal = 'TRBL8';

singleAnimal = false;               % processing of single animal or average ?

% animalTimeWidth   = 20;
% animalTimeOverlap = 0;  % to make a plot at every single frame
animalTimeWidth   = 2;
animalTimeOverlap = 0.96;  % to make a plot at every single frame

% multiTimeStart    = '00h00';
% multiTimeStop     = '20h00';
multiTimeStart    = '12h00';
multiTimeStop     = '32h00';    % 28 for wt3

%%% Time related (mod 2.1)
temperature = 25;       % to correct "delta_t" (2.1)
delta_t = 5;                        % LTA (moved here 2.1)
timeShift = 0;                  % in HOURS, correction of reference time (2.1)
% NB: will be stored in TIME structure

%%% Grid related:
gridType = 'L';         % 2.3
gridOverlap = 0.5;      % 2.3
% TO PROCESS JESUS MOVIES (26h00 olap 0)
% gridTime = '26h00';               % time APF used to draw grid and determine cell patches to track (2.1)
% gridTimeLandmarks = '26h00';      % SHOULD DISAPPEAR AS ALL GRIDS SHOULD BE DRAWNED AT SAME TIME AND THERE SHOULD ONLY BE "gridTime"
% TO PROCESS BIG WT MOVIES (128 pix, olap 0.5)
gridTime = 'start';                 % time APF used to draw grid and determine cell patches to track (2.1)
gridTimeLandmarks = '16h30';        % SHOULD DISAPPEAR AS ALL GRIDS SHOULD BE DRAWNS AT SAME TIME AND THERE SHOULD ONLY BE "gridTime"


% Quantities to be PLOTTED and PROJECTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qname = {'EG';'S'};  % Quantity names 'Rho','M','CDcad','EG',dnA...
Qname = {'EG';'dnA';'S'};  % Quantity names 'Rho','M','CDcad','EG',dnA...
% Qname = {'EG';'dnA';'S'};  % Quantity names 'Rho','M','CDcad','EG',dnA...
% Qname = {'EA';'ER';'ES';'ED';'EG';'Phi';'dnA';'dnD';'S'};  % Quantity names 'Rho','M','CDcad','EG',dnA...
% Qname = {'EA';'ER';'ES';'ED';'EG';'S';'M';'I'};  % Quantity names
% Qname = {'dnA';'dnD'}; 
% NB: ALL quantities are CALCULATED anyway in AOA, DOA, *BUT* in POA, ONLY those listed here will be calculated and plotted in POA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Projection (moved 2.1) %%

% NB: relevant for POA, LTA

% if projection is involved, specify ONTO what here below (1st run POA, then LTA):
% uAnimal = 'meanTRBL_grid26h';
% uAnimal = 'TRBL8';
% uAnimal = 'R23';
% uAnimal ='meanWT_grid26h';
uAnimal ='meanWT';
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
PLOT.boxSize = 128*[1 1];           % box size of the average grid (BIG and TRBL movies)
% PLOT.boxSize = [122 128/2]/0.322;   % JESUS origin & grid size FOR HALF MOVIES
% PLOT.boxSize = [122 128]/0.322;   % JESUS origin & grid size FOR FULL MOVIES
PLOT.macrocaetes = 1;               % plot macrocaetes position
PLOT.displayOrigin = true;          % 2.6
PLOT.drawMidline = true;            % 2.3
PLOT.makeItFull = true;             % only applies when processing HALF notums (2.3)

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
allAnimals = { 'BIG_1' ; 'BIG_5' ; 'BIG_6' ; 'BIG_N5m' ; 'BIGwt2r' ; 'BIGwt2l'; 'TRBL1' ; 'TRBL4' ; 'TRBL7' ; 'TRBL9' ; 'TRBL8l'  ; 'TRBL8r'  };
% allAnimals = { 'TRBL1' ; 'TRBL4' ; 'TRBL7' ; 'TRBL9' ; 'TRBL8l'  ; 'TRBL8r'  }; %#ok<*UNRCH> % animal name list
% allAnimals = { 'BIG_1' ; 'BIG_5' ; 'BIG_6' ; 'BIG_N5m' ; 'BIGwt2r' ; 'BIGwt2l'};
% allAnimals = {'BIGwt2'};
% allAnimals = {'wt3'};
%------------------------------------------------------------

% Animals to take into account in the average:
%------------------------------------------------------------
% avgAnimals = { 'TRBL1' ;'TRBL4' ; 'TRBL7' ; 'TRBL9' ; 'TRBL8l' ; 'TRBL8r'};
avgAnimals = { 'BIG_1' ; 'BIG_5' ; 'BIG_6' ; 'BIG_N5m' ; 'BIGwt2r' ; 'BIGwt2l'};
% avgAnimals = {'BIGwt2'};
% avgAnimals = {'wt3'};

% ESCARGOT FULL 12 APF
% allAnimals = {'esgGfp_12apf_1';'esgGfp_12apf_2';'esgGfp_12apf_3';'esgGfp_12apf_4'};
% avgAnimals = {'esgGfp_12apf_1';'esgGfp_12apf_2';'esgGfp_12apf_3';'esgGfp_12apf_4'};

% ESCARGOT FULL 24 APF
% allAnimals = {'esgGfp_24apf_1';'esgGfp_24apf_2';'esgGfp_24apf_3'};
% avgAnimals = {'esgGfp_24apf_1';'esgGfp_24apf_2';'esgGfp_24apf_3'};
%------------------------------------------------------------

% PnameAOA = 'TA'; % Process name: SM,VM,AOS,TA ** TO BE REMOVED **
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

% type of projection to be PLOTTED in maps (only one at a time)
tagPlot = 'd';  % choose between 'd', 'i', 'do' for, respectively deviator, isotropic, orthogonal deviator

% project this:
OnamePOA = 'AOA';              % Operation name (AOA,DOA...) (overridden by "empty" when single animal processing)

QunitsPOA = 'h^{-1}';
QsrPOA    = 50*10^1;
QsbPOA    = 0.1;



%% LTA - Local Time Analysis %%

% Enter AOA to average mean quantities Q such as Rho,dnA... If Q is a tensors (EG,...) will plot |Q| (tensor norm)
% Enter POA and "QdotuQ" to plot the projection of Q along uQ
% NB: units used in the axes are the ones corresponding to LAST item listed in "QnameLTA"

OnameLTA = 'POA';         % Operation name (AOA, DOA, POA...)
 

QnameLTA = {'SdotEG','EGdotEG'}; % list of quantities to be PLOTTED
Qfactor =  [ 7.5e-1      1]; % to be able to plot different type of data (2.3)

% QnameLTA = {'dnA','SdotEG','EGdotEG'}; % list of quantities to be PLOTTED
% Qfactor =  [  0.1   7.5e-1      1]; % to be able to plot different type of data (2.3)

% QnameLTA = {'ERdotEG','ESdotEG','EDdotEG','EAdotEG','SdotEG','MdotEG','IdotEG','EGdotEG'}; % list of quantities to be PLOTTED
% Qfactor =  [    1        1           1       1       7.5e-1    3.2e-3    6.4e-2     1]; % to be able to plot different type of data (2.3)
%     QnameLTA = {'EG','ED','ER','ES','EA','M','I','Rho'};
%     QnameLTA = {'EDdotEG','ERdotEG','ESdotEG','EAdotEG'};
%     QnameLTA =  {'dnA','dnD'};
%     Qfactor =   [  1   , 1e-1]; % to be able to plot different type of data (2.3)

% range of zone xy coordinates (not ij!!!) to be included (those coordinates start at [1 1]):
Qareas = {[5 13; 1 7]}; % {[xMin xMax ; yMin yMax]} of box coodinates         % meanWT
%     Qareas = {[2 2; 2 2]}; % {[xMin xMax ; yMin yMax]} of box coodinates      % TRBL8
% Qareas = {[2 2; 1 1]}; % {[xMin xMax ; yMin yMax]} of box coodinates            % 26h00: meanTRBL & meanWT
%     Qareas = {[1 1; 1 1]}; % {[xMin xMax ; yMin yMax]} of box coodinates      % R23
%     Qareas = {[2 2; 3 3]}; % {[xMin xMax ; yMin yMax]} of box coodinates      % BIGwt2
% NB: Examples: {[1 6; 1 5]} will process all boxes from (1,1) to box (6,5) (in xy coordinates)
% NB: Examples: {[6 6; 5 5]} will ONLY process the box (6,5) (in xy coordinates)

% Axis range, leave empty for automatic range selection, otherwise like this: {[400 410]}
    Qrange = {};
% Qrange = {0.07*[-1 1]}; % BIG & TRBL
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

nQname = length(Qname); % 2.2

% defining fullTag (2.3)
fullTag = '';
if PLOT.makeItFull
    fullTag = '_full';
end

% Format, resolution, and extention
printFormat     = ['-d' PLOT.extension];
printResolution = ['-r' num2str(PLOT.resolution)];
imageExtension  = ['.'  PLOT.extension];

% Store time information in struct more simple managment
TIME.animalTimeOverlap = animalTimeOverlap;
TIME.animalTimeWidth   = animalTimeWidth;
TIME.multiTimeStart    = multiTimeStart;
TIME.multiTimeStop     = multiTimeStop;
TIME.gridTime = gridTime;       % 2.2
TIME.delta_t = delta_t;         % 2.2
TIME.gridType = gridType;       % 2.3
TIME.gridOverlap = gridOverlap; % 2.3

CustomColors;       % defines usual set of colors
AllQsColorsUnits;   % associate quantities with specific colors

% display parameters to be used during plot
minimalInfoDisplay = false;
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

% Start filling DISPLAY (2.3)
DISPLAY.makeItFull = PLOT.makeItFull;
DISPLAY.displayOrigin = PLOT.displayOrigin; 
DISPLAY.plotType = QplotType;
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
DISPLAY.fadeColor = custom_white; % dont know why but plotfield need it for VM
DISPLAY.drawMidline = PLOT.drawMidline;
DISPLAY.midlineColor = grey;

% Scale values according to the information to be ploted
sr_AOS =         {8e3 ;  5e2  ; 5e2 ; 25 ;  1 ;  5  ; 50 ; 4e3  ;  8e3   ;  [10 70]};   % sets ratio setting size of ellipses or bars in tensor representation for M,I.
srbar_AOS =      [1e-2 ; 1e-1 ; 0.1 ;  2 ; 50 ; 10  ;  1 ; 5e-2 ;  1e-2  ;   2];   % scale bar lengths for each contribution     
killtrace_AOS =  [ 0   ;  0   ;  0  ;  0 ;  1 ;  0  ;  0 ;  0   ;    0   ;    1];
% allQs_AOS = BoxCellArea Rho    I     M    V   dnD  dnA rRho rBoxCellArea CDcad CDesg CDmyo CDsqh % NB: SHOULD ALWAYS BE LISTED IN THIS ORDER!!

% sr_AOS        = { 4e2 ; 50 ;0.8*[1 3]; 100 ; [4 30] };   % sets ratio setting size of ellipses or bars in tensor representation for M,I. % 2 hours
% % sr_AOS        = { 5e2 ; 50 ;[1 3]; 100 ; [4 30] };   % sets ratio setting size of ellipses or bars in tensor representation for M,I. % 14h/20h hours
% srbar_AOS     = [ 0.1 ;  2 ; 50 ;  1  ;    10  ];   % scale bar lengths for each contribution
% killtrace_AOS = [  0  ;  0 ;  0 ;  1  ;     1  ];   % Will set average compartment trace to 0 in the plots (mean isotropic part = 0). Choose this when tensors are known up to an additive constant
% % contributions   Rho    I    M    V       CD

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
sr_TA        = {  5e3 ;  2e3 ;  5e3 ;  8e3 ;  2e3 ;  9e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  4e3 ;  6e3 ;  4e3 };  % ISO average 14h/20h
% sr_TA        = 5e3; sr_TA = num2cell(sr_TA*ones(15,1));  % DEV average 14h
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
% contributions     G      S      R     Ds      D      A      N      F      J     Jb     DM      E    EPSI    Phi   Egeo
 
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
    disp('Load AIA_info of each animal and determine common AreaRatios ...')
    
    % Values that will NOT be averaged
    rejectListOld = {'xywh';'size';'overlap';'centroids';'ULCs'; 'fullImage';'color';'lineWidth';'coordinates'; 'TimeArray';'FrameArray'};
    rejectListNew = {'EDM'; 'EGeig'; 'EPSIeig'; 'ESeig'; 'Eeig' ; 'Phieig' ;'RCondsRaw'}; % 2.2
    rejectListAreaRatios = {'AreaRatios_TA';'AreaRatios_SM';'AreaRatios_AOS';'AreaRatios_GEP'; 'RConds'}; % will only average common "AreaRatios" (2.2)
    rejectList = [rejectListOld ; rejectListNew ; rejectListAreaRatios]; % 2.2
    
    nAllAnimals = length(allAnimals);       % 2.2
    nAvgAnimals = length(avgAnimals);       % 2.2
    backupPathList = {};                    % initialisation for list of path to the mean backup
    macroCoord = NaN(16,2,nAllAnimals); % upped from 14 (2.1)
    deltaYmid = NaN(nAllAnimals,1);     % will store RESCALED distance between midline and origin IN PIXELS (2.3)
    
    % loading parameters and data for each animal using their respective AIA_info_
    mean_GRID_specs = cell(nAllAnimals, 1);
    for n = 1:nAllAnimals
        
        disp(allAnimals{n});
        
        [GRID backupPath RAW] = SingleAnimalLoader(AIAFolderName, allAnimals{n}, TIME, true );  % removed PnameAOA (2.2), moved out of if (2.3)
        
        % Checking animal is not full already (2.3)
        if PLOT.makeItFull && strcmp(RAW.halfNotum,'b')
            warndlg('Cannot have "makeItFull=true" when processing a full animal!! Please set it to "false".','ERROR!');
            return
        end
        
        %%% determines RESCALED deltaYmid for this animal (2.3)
        if ~isempty(RAW.yMid)
            deltaYmid(n) = (RAW.xyStart(2) - RAW.yMid)*RAW.yFactor; % directly applying factor
        end
        
        % If it is an animal to be processed in the mean
        if ismember( allAnimals{n}, avgAnimals )
            [~,cpt] = ismember(allAnimals{n}, avgAnimals); % for the AOS multichannel            
%             [GRID backupPath AOSvariable, RAW] = SingleAnimalLoader(AIAFolderName, allAnimals{n}, PnameAOA, TIME, true );
            
            % Store backup path to be loaded later
            backupPathList = [backupPathList backupPath]; 
            thisAOTbackup = load(backupPath);

            %%% Determining common AreaRatios giving priority to AreaRatios based on segmented (images)(2.2)
            if isfield(thisAOTbackup,'AreaRatios_TA')
                AreaRatios = thisAOTbackup.AreaRatios_TA;
                RcondsNoNaN = thisAOTbackup.RConds;       % initializes RcondsNoNaN with RConds
                RcondsNoNaN(isnan(RcondsNoNaN)) = 0;      % Replacing NaNs by 0 because AreaRatios must have ONLY 0s and no NaNs (where no animal)
                AreaRatios = AreaRatios.*RcondsNoNaN;     % including RConds in TA AreaRatios
            elseif isfield(thisAOTbackup,'AreaRatios_SM') % SM before AOS because SM directly eliminates border cells
                AreaRatios = thisAOTbackup.AreaRatios_SM;
            elseif isfield(thisAOTbackup,'AreaRatios_AOS') 
                AreaRatios = thisAOTbackup.AreaRatios_AOS;
            elseif isfield(thisAOTbackup,'AreaRatios_GEP') % last because also exists when not segmented
                AreaRatios = thisAOTbackup.AreaRatios_GEP;
            end
            
            % creates "GRID_AOT_animal" structure for each animal (2.2)
            thisAOTbackup.AreaRatios = AreaRatios;
            eval(['GRID_AOT_' allAnimals{n} ' = thisAOTbackup;']);% 2.2
            
            %%% loads macrocaete position
            %--------------------------------------------------------------
            nbMacro = 8;
            if strcmp(RAW.halfNotum,'b')
                nbMacro = 16; % upped from 14 (2.1)
            end
            macroCoord(1:nbMacro,:,n) = LandmarkLoader(RAW.rawPathFolder, allAnimals{n}, RAW.halfNotum, gridTimeLandmarks); % added "gridTime" (2.1), replaced by "gridTimeLandmarks" (2.3)
%             macroCoord(1:nbMacro,:,n) = LandmarkLoader(RAW.rawPathFolder, allAnimals{n}, RAW.halfNotum);
            coefXY = PLOT.boxSize ./ RAW.boxSize;
            macroCoord(:,1,n) = (macroCoord(:,1,n) - RAW.xyStart(1)) .* coefXY(1);
            macroCoord(:,2,n) = (macroCoord(:,2,n) - RAW.xyStart(2)) .* coefXY(2);
            %--------------------------------------------------------------
            
            % COMMENTED 2.3
%         else
%             [GRID, ~, ~] = SingleAnimalLoader(AIAFolderName, allAnimals{n}, TIME, false ); % removed PnameAOA (2.2)
        end
        
        % save the grid in the mean grid structure
        mean_GRID_specs{n} = GRID;
        gridOverlap = GRID.overlap;
    end

    macrocaetes = NaN(5,16);
    macrocaetes(1:2,:) = nanmean(macroCoord,3)';
    DISPLAY.macrocaetes = macrocaetes;
    meanDeltaYmid = nanmean(deltaYmid); % 2.2
    
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
        %         [GRID backupPathList] = SingleAnimalLoader(AIAFolderName, deltaAnimals{1}, PnameDOA, TIME, true );
        %         disp(backupPathList)
        %         load(backupPathList);
        %         LoadedBackup{1} = eval(['GRID_' PnameDOA '_' deltaAnimals{1}]);
        %
        %         % animal 2
        %         [GRID backupPathList] = SingleAnimalLoader(AIAFolderName, deltaAnimals{2}, PnameDOA, TIME, true );
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
        [~, backupPathList, RAW] = SingleAnimalLoader(AIAFolderName, Animal,  'AOT',  TIME,  true );
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
    
    tagTensor = {'i';'d';'do'}; %#ok<*UNRCH>
    tagProjectionTime = [strrep(num2str(uTimeWidth),'.','') '_' strrep(num2str(uTimeOverlap),'.','')];
    
    QplotType = 'circle';
    QKillTr = 0;
    
%     Pname  = GetPname(Qname{1});
%     uPname = GetPname(uQname  );
    
    uTIME = TIME;
    uTIME.animalTimeWidth = uTimeWidth;
    uTIME.animalTimeOverlap = uTimeOverlap;
    
    if isempty(OnamePOA)
        
        % projection is a SINGLE animal: 
        % - Animal must be the same as uAnimal
        % - OnamePOA and uOname must be empty
        % - Qname can be different than uQname
        
        uOname = '';
        uAnimal = Animal;
        PLOT.significance = false;
        
        [~, backupPathList, RAW] = SingleAnimalLoader(AIAFolderName, Animal, TIME, true ); % 2.3
        [~, ubackupPathList]  = SingleAnimalLoader(AIAFolderName, uAnimal, uTIME, true );  % 2.3
%         [~, backupPathList, RAW] = SingleAnimalLoader(AIAFolderName, Animal,  'AOT',  TIME, true );
%         [~, ubackupPathList]  = SingleAnimalLoader(AIAFolderName, uAnimal, 'AOT', uTIME, true );
        
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
        
        % projection is an AVERAGED animal: 
        % - check if uAnimal is specified, else uAnimal = Animal
        % - check if uOname is specified, else uOname  = OnamePOA
        % - Qname can be different than uQname
        
        if isempty(uAnimal), uAnimal = Animal; end
        if isempty(uOname),  uOname  = OnamePOA;  end
        
        animalProjectionFolder = [OnamePOA  '_' Animal  '_' num2str(TIME.animalTimeWidth)  'h_olap_' num2str(TIME.animalTimeOverlap) ];
        animalProjectorFolder  = [uOname '_' uAnimal '_' num2str(uTIME.animalTimeWidth) 'h_olap_' num2str(uTIME.animalTimeOverlap)];
        
        disp( 'Loading backup TO BE projected...' )
        thisFilename = [PathName filesep animalProjectionFolder filesep OnamePOA '_backup'];
        disp( thisFilename );
        BACKUP = load(thisFilename);
        
        disp( ['Loading backup FOR projection...' uQname] )
        thisFilename = [PathName filesep animalProjectorFolder filesep OnamePOA '_backup'];
        disp( thisFilename );
        TEMPBACKUP = load( thisFilename );
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
            [~, backupPathList, RAW] = SingleAnimalLoader(AIAFolderName, Animal,  'AOT',  TIME,  true );
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

% 18-26/01/2017: 2.3 (Boris)
% - added calculation of RESCALED mean y distance between origin and midline: "meanDeltaYmid"
% - changes to adapt to latest versions of AOA (2.7), POA (3.3) and LTA(3.4)
% - added "Qfactor" in LTA to be able to scale separately each quantity

% 12,18/01/2017: 2.2 (Boris)
% - ONLY using AOT backups now!! (=> removed "PnameAOA"...)
% - removed "SingleAnimalLoader" output "AOSvariable"

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
