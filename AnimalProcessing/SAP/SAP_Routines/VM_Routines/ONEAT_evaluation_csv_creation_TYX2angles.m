
%% ONEAT_evaluation_csv_creation_TYX2angles.m %%


% Creation of CSV files for ONEAT to evaluate MATLAB predictions (Apoptose and Divisions (with angles).
% Order of coordinates : T Y X JunctionAngle CentroïdAngle
% Lucas Sancéré




%% Select Animal %%

animalName = 'BIG1';

%% Load informations in SAP_info and CPT

MAPcall = true;
eval(['SAP_info_' animalName ';']);
number_of_frame = finalFrame;

load([pathFolderCTD filesep 'allDelaminatingCells.mat']);
load([pathFolderCTD filesep 'allDividingCells.mat']);
load([pathFolderCTA filesep 'CTA_' animalName '.mat']) 

%% Creation of csv for bridge

%Table.mat for delaminations

Coordinates_Delaminations = round(allDelaminatingLastXYs);
T_Y_X_table_Del = zeros(length(allDelaminatingLastXYs),3);
T_Y_X_table_Del(:,1) = allLastFramesDel;  % allLastFramesDel comes from allDelaminatingCells.mat from CTD
T_Y_X_table_Del(:,2) =  Coordinates_Delaminations(:,2);
T_Y_X_table_Del(:,3) =  Coordinates_Delaminations(:,1);


%Table.mat for divisions

Coordinates_Divisions = round(allDividingLastXYs);
T_Y_X_alpha1_alpha2_table_Div =  zeros(length(allDividingLastXYs),5);
T_Y_X_alpha1_alpha2_table_Div(:,1) = allLastFramesDiv;  %  allLastFramesDiv comes from allDividingCells.mat from CTD
T_Y_X_alpha1_alpha2_table_Div(:,2) = Coordinates_Divisions(:,2);
T_Y_X_alpha1_alpha2_table_Div(:,3) = Coordinates_Divisions(:,1);
T_Y_X_alpha1_alpha2_table_Div(:,4) = allJunctionAngles(:,1);
T_Y_X_alpha1_alpha2_table_Div(:,5) = allCentroidAngles(:,1);




%% Conversion of files and saving

% Creation of Folders

 pathCPT_ONEAT_Input = [gridFolderCPT filesep 'ONEAT_Input'];
 mkdir(pathCPT_ONEAT_Input)
 
 gridFolderONEAT_Input = [pathCPT_ONEAT_Input filesep];
 pathCPT_ONEAT_Input_Del = [gridFolderONEAT_Input filesep 'Delaminations_Coordinates_TYX'];
 mkdir(pathCPT_ONEAT_Input_Del)
 pathCPT_ONEAT_Input_Div = [gridFolderONEAT_Input filesep 'Divisions_Coordinates_TYXalphas'];
 mkdir(pathCPT_ONEAT_Input_Div)

 
% Save files into this folders

 save([pathCPT_ONEAT_Input_Del filesep animalName '_Delaminations_Coordinates_TYX.mat'], 'T_Y_X_table_Del');
 save([pathCPT_ONEAT_Input_Div filesep animalName '_Divisions_Coordinates_TYXalphas.mat'], 'T_Y_X_alpha1_alpha2_table_Div');
 

% Convert into csv files 

 T_Y_X_csv_Del = load([pathCPT_ONEAT_Input_Del filesep animalName '_Delaminations_Coordinates_TYX.mat']);
 csvwrite( 'Delaminations_Coordinates_TYX.csv', T_Y_X_csv_Del.T_Y_X_table_Del);
 
 T_Y_X_alpha1_alpha2_csv_Div = load([pathCPT_ONEAT_Input_Div filesep animalName '_Divisions_Coordinates_TYXalphas.mat']);
 csvwrite( 'Divisions_Coordinates_TYXangles.csv', T_Y_X_alpha1_alpha2_csv_Div.T_Y_X_alpha1_alpha2_table_Div);

%% History

%10/11/2020: 1.0 (Lucas)
% -Creation
