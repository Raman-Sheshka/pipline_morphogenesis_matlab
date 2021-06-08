% Loads two AOA backup and make the difference between the two
%
% Stephane Rigaud
% version 1.3
%


%% DOA Process

% list information relative to the grid, not concerned by the operation, to be copied without modification
gridInfo = {'xywh';'size';'overlap';'color';'lineWidth';'fullImage';'centroids';'ULCs';'TimeArray';'FrameArray';'coordinates'};
% list information concerned by the operation but to be treated differently
exceptions = {'RConds';'Macrocaetes';'AreaRatios';'AreaRatios_TA';'AreaRatios_SM';'AreaRatios_VM';'AreaRatios_AOS';'errorDnPs';'errorPs'};
% all other information contained in the backup will be processed as follow: A.x - B.x = C.x

DOA_backup = struct(); % prepare the new backup structure
fieldList = fieldnames(LoadedBackup{1}); % Get the list of field to process from one of the backup
for i = 1:length(fieldList)
    if isfield(LoadedBackup{2},fieldList(i))
        
        if ismember(fieldList(i), gridInfo) % if grid information, simply copy of the field
            disp(['grid info : ' fieldList{i}]);
            eval(['DOA_backup.' fieldList{i} ' = LoadedBackup{1}.' fieldList{i} ';']);
        else
            if ismember(fieldList(i),exceptions) % if exception value, treat the value accordingly
                if strcmp(fieldList{i},'Macrocaetes') % mean(A.macro, B.macro)
                    disp(['mean : ' fieldList{i}]);
                    DOA_backup.Macrocaetes = (LoadedBackup{1}.Macrocaetes + LoadedBackup{2}.Macrocaetes) ./ 2;
                elseif strcmp(fieldList{i},'RConds') % min(A.RConds, B.RConds)
                    disp(['min : ' fieldList{i}]);
                    DOA_backup.RConds = nanmin(LoadedBackup{1}.RConds, LoadedBackup{2}.RConds);
                elseif strncmp(fieldList{i},'AreaRatios_',11) % min(A.ar, B.ar)
                    disp(['min : ' fieldList{i}]);
                    eval(['DOA_backup.' fieldList{i} ' = min(LoadedBackup{1}.' fieldList{i} ',LoadedBackup{2}.' fieldList{i} ');']);
                elseif strncmp(fieldList{i},'error',5) % max(A.error, B.error)
                    disp(['max : ' fieldList{i}]);
                    eval(['DOA_backup.' fieldList{i} ' = max(LoadedBackup{1}.' fieldList{i} ',LoadedBackup{2}.' fieldList{i} ');']);
                end
            else
                if ~isempty(strfind(fieldList{i},'_Smap')) % if Significance map, logical AND(A.map, B.map)
                    disp(['AND : ' fieldList{i}]);
                    eval( ['DOA_backup.' fieldList{i} ' = or( LoadedBackup{1}.' fieldList{i} ', LoadedBackup{2}.' fieldList{i} ');'] );
                else
                    disp(['Substract : ' fieldList{i}]);
                    eval( ['DOA_backup.' fieldList{i} ' = LoadedBackup{1}.' fieldList{i} ' - LoadedBackup{2}.' fieldList{i} ';'] );
                end
            end
        end
        
    end
end

%%% save results
disp('Save DOA backup ...')
if ~exist(OutputFolderName,'dir'); mkdir(OutputFolderName); end
save( [OutputFolderName filesep 'DOA_backup'], ['DOA_backup']);




%% plot
if ~isempty(Qname) && PLOT.plot
    
    %%% reload the backup (not necessary, we still do it as a verification process)
    load([OutputFolderName filesep 'DOA_backup']);
    GRID = DOA_backup;
    
    %%% check dimension (compatibility for single time point plot)
    Lt = size( eval(['GRID.' Qname{1}]) ,4);
    if Lt == 0
        Lt = size( eval(['GRID.' Qname{1}]) ,3);
    end
    
    %%% initialisation of the plot
    % get the BG image (simple white image)
    size_image_x = PLOT.boxSize(1) * (GRID.size(2) * GRID.overlap + 1);
    size_image_y = PLOT.boxSize(2) * (GRID.size(1) * GRID.overlap + 1);
    if GRID.overlap == 0
        size_image_x = PLOT.boxSize(1) * (GRID.size(2) + 1);
        size_image_y = PLOT.boxSize(2) * (GRID.size(1) + 1);
    end
    Qimage = ones([size_image_y , size_image_x]);
    
    %%% get the display info parameters
    DISPLAY.Animal             = DOAname;
    DISPLAY.plotType           = QplotType;
    DISPLAY.minimalInfoDisplay = minimalInfoDisplay;
%     DISPLAY.minAEV             = minAEV;
    DISPLAY.scaleBarWidth      = scaleBarWidth;
    DISPLAY.gridDisplay        = gridDisplay;
    DISPLAY.lineWidth          = lineWidth;
    DISPLAY.pointSize          = pointSize;
    DISPLAY.EVstyles           = EVstyles;
    DISPLAY.signOpacities      = signOpacities;
    DISPLAY.fontSizeInfo       = fontSize;
    DISPLAY.imageFading        = imageFading;
    DISPLAY.errorPsMin         = 1; % 10^-10; % hard coded, todo
    DISPLAY.errorDnPsMin       = 1; % 10^-10;
    DISPLAY.errorFontSize      = 5;
    if PLOT.macrocaetes
        DISPLAY.macrocaetes    = GRID.Macrocaetes;
    end
    DISPLAY.fadeColor = custom_white; % dont know why but plotfield need it for VM


    for t = 1:Lt
        for j=1:length(Qname)
            
            if PLOT.Qparameters
                [~,idx] = GetPname(Qname{j});
                if idx > 0
                    Qcolor  = eval(['allColors_' Pname '{idx}']);
                    Qunits  = eval(['allUnits_' Pname '{idx}']);
                    Qsr     = eval(['sr_' Pname '{idx}']);
                    Qsrbar  = eval(['srbar_' Pname '(idx)']);
                    QKillTr = eval(['killtrace_' Pname '(idx)']);
                end
            end
            eval(['GRID.AreaRatios = GRID.AreaRatios_' Pname ';']);
            
            % kill trace tag
            KillTrTag = '';
            if QKillTr
                KillTrTag = '_Tr=0';
            end
            
            % Set the time and n variable
            DISPLAY.time = [ GRID.TimeArray{t,1} ' - ' GRID.TimeArray{t,2} ];
            DISPLAY.n = round( ( GRID.FrameArray(t,1) + GRID.FrameArray(t,2) ) / 2 );
            
            if Lt > 1 %(1.2)
                DISPLAY.step = t;
            end
            
            % significant_map_wt = Significativity(D,dT_wt,t);
            if PLOT.significance
%                 if PLOT.error % in progress v2
%                     stdT = eval( [Oname '_stdMAP.' Qname{j} '_std;'] );
%                     meanT = eval( ['GRID.' Qname{j} ';'] );
%                     temp_sign_map = ErrorSignificativityMap( meanT, stdT);
%                 else
                    temp_sign_map = eval( ['GRID.' Qname{j} '_Smap;'] );
%                 end
                DISPLAY.significant_map = temp_sign_map(:,:,:,t);
            end
            
            % Plot image
            PlotField(Qname{j},QKillTr,Qcolor,Qunits,Qsr,Qsrbar,GRID,Qimage,DISPLAY);
            
            % Print plot
            if PLOT.print
                this_filename = ['DOA_' DOAname '_' Qname{j} '_' GRID.TimeArray{t,1} 'to' GRID.TimeArray{t,2} KillTrTag '_sr=' num2str(Qsr(1)) printFormat];
                if ~exist([OutputFolderName filesep 'Frames_DOA_' QplotType],'dir'); mkdir([OutputFolderName filesep 'Frames_DOA_' QplotType]); end
                this_filename = ['DOA_' DOAname '_' Qname{j} '_' GRID.TimeArray{t,1} 'to' GRID.TimeArray{t,2} KillTrTag '_sr=' num2str(Qsr(1)) printFormat];
                if strcmp(imageExtension,'.svg')
                    plot2svg([OutputFolderName filesep 'Frames_DOA_' QplotType filesep this_filename],figure(1),'png');
                else
                    print (printFormat, printResolution, [OutputFolderName filesep 'Frames_DOA_' QplotType filesep this_filename]);
                end
            end
            close
        end % end for each quantities
    end % end for each time point
end



%% historique
 
% 19/03/2015: 1.3
% - code style modification and improvement
% - manage the significativity in boolean map

% 21/01/2015: 1.2
% - fixed bug when nb time point data == 1

% 13-01-2015 : 1.1
% - replot and print option
% - make difference between all the values in backups



