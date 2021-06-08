function [GRID backupPathList AOSvariable RAW] = SingleAnimalLoader(AIAFolderName, Animal, Pname, TIME, tobeprocessed )

AIA_call = 1;

% return
backupPathList = [];
AOSvariable = struct();

% run AIA_info to retrive ... info
run([ AIAFolderName filesep 'AIA_info_' Animal ]);

% rebuild the GRID
GRID = GridMaker(imageSize, boxSize, xyStart, gridSize, gridColor, gridLineWidth, gridOverlap); % added "gridOverlap"
xywh = GRID.xywh;
nx   = GRID.size(2);
ny   = GRID.size(1);
gridOverlap = GRID.overlap; % has been set to 0 by "GridMaker" if was empty in AIA_info
olapTag = ''; % default value: not displaying overlap in folder name when empty or 0 (for compatibility)
if gridOverlap > 0
    olapTag =  ['_olap_' num2str(gridOverlap)];
end

RAW.nDigits = nDigits;
RAW.filename = filename;
RAW.rawPathFolder = pathFolderRaw;
RAW.gridAnimalFolder = '';
RAW.imageFormat = imageFormat;
RAW.boxSize = boxSize;
RAW.xyStart = xyStart;
RAW.halfNotum = halfNotum;


if tobeprocessed
    
    % Defines grid specific subfolder name
    gridSpecsName = ['Grid_xy_' num2str(xywh(1)) '_' num2str(xywh(2)) '_wh_' num2str(xywh(3)) 'x' num2str(xywh(4)) '_nynx_' num2str(ny) 'x' num2str(nx) olapTag];
    
    %%% Defines the "Process_Folder", "Grid_folder", "Average folder"
    %--------------------------------------------------------------
    % Process
    backupPath = eval(['pathFolder_' Pname]);
    %--------------------------------------------------------------
    
    %--------------------------------------------------------------
    % Grid
    tag = '';
    if strcmp(Pname,'SM')
        tag = ['_' code2run '_dmu=' num2str(muAccuracy)]; % specific case for SM (for compatibility)
    end
    backupPath = [backupPath filesep gridSpecsName];
    gridLBackupFolder = [ pathFolder_TA filesep gridSpecsName ];  % path to the grid backups containing the Lgrid information
    Grid_backup_rootfilename = eval(['filename_' Pname]);      % 1.9
    RAW.gridAnimalFolder = gridLBackupFolder;
    %--------------------------------------------------------------
    
    %--------------------------------------------------------------
    % Average
    % if AOA_time start and stop are define use them, else use time defined in AOT of AIA_parameters
    averageFolderName = ['Average_' num2str(TIME.animalTimeWidth) 'h_' TIME.multiTimeStart '_to_' TIME.multiTimeStop '_olap_' num2str(TIME.animalTimeOverlap)]; % moved down here (1.9)
    backupPath = [backupPath filesep averageFolderName];
    %--------------------------------------------------------------
    
    
    %%% Loads "mean" backup
    % NB: multiple name template possible depending on the backup:
    % mean_PROCESS_ANIMAL.mat | mean_PROCESS_ANIMAL_Backup.mat | mean_GRIDTYPE_PROCESS_ANIMAL.mat
%     fullpath_backup_PIV = [ backupPath filesep 'mean_' Pname '_PIV_' Animal '_Backup.mat'];
    fullpath_backup     = [ backupPath filesep 'mean_' Pname '_' Animal '_Backup.mat'];
    fullpath_lagrangian = [ backupPath filesep 'mean_' Pname '_LGrid_' Animal '.mat'];
    fullpath_eulerian   = [ backupPath filesep 'mean_' Pname '_EGrid_' Animal '.mat'];
    fullpath            = [ backupPath filesep 'mean_' Pname '_' Animal '.mat'];
    if exist(fullpath_backup,'file')
        load(fullpath_backup);
        backupPathList = [fullpath_backup];
    elseif exist(fullpath,'file')
        load(fullpath);
        backupPathList = [fullpath];
    elseif exist(fullpath_lagrangian,'file')
        load(fullpath_lagrangian);
        backupPathList = [fullpath_lagrangian];
    elseif exist(fullpath_eulerian,'file')
        load(fullpath_eulerian);
        backupPathList = [fullpath_eulerian];
%     elseif exist(fullpath_backup_PIV,'file')
%         backup = load(fullpath_backup_PIV);
%         backupPathList = [fullpath_backup_PIV];
    else
        disp(['Error: backup "' fullpath '" was not found and was skipped.'])
        return;
    end
        
    %%% Store AOS multichannel name, to be used later (compatibility)
    AOS_sub_name = [];
    AOS_dif_name = [];
    AOS_tmp_name = [];
    if strcmp(Pname,'AOS')
        if length( filenameRaw ) > 1
            idx = find( ismember( filenameRaw{1}, '_' ) );
            for r = 1:length( filenameRaw )  % for each signal 1 a n
                AOS_tmp_name{r} = filenameRaw{r}(1:idx-1); % get the signal name eg {'cad1_' ; 'sqh1_'} => cad1 or sqh1
            end
            AOS_dif_name = intersect( AOS_tmp_name{1}, AOS_tmp_name{end}, 'stable' );
        else
            AOS_sub_name = filenameRaw{1};
        end
    end
    
    AOSvariable.subName = AOS_sub_name;
    AOSvariable.difName = AOS_dif_name;
    AOSvariable.tmpName = AOS_tmp_name;
    
end

end