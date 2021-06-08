% DifferenceOverAnimals (DOA)
% 
% Loads two AOA backup and make the difference between the two
%
% version 1.5
% Stephane Rigaud
% Boris Guirao



%% DOA Process

% list information relative to the grid, not concerned by the operation, to be copied without modification
% gridInfoList = {'xywh';'size';'overlap';'color';'lineWidth';'fullImage';'centroids';'ULCs';'TimeArray';'FrameArray';'coordinates'};
% list information concerned by the operation but to be treated differently
% exceptions = {'RConds';'Macrocaetes';'AreaRatios';'AreaRatios_TA';'AreaRatios_SM';'AreaRatios_VM';'AreaRatios_AOS';'errorDnPs';'errorPs'};
% all other information contained in the backup will be processed as follow: A.x - B.x = C.x

gridInfoListExt = [gridInfoList ; 'AreaRatios'; 'errorDnPs'; 'errorPs']; % adding those in gridInfoList defined in AIA_MultiOperation (3.6)

DOA_backup = struct(); % prepare the new backup structure
fieldList = fieldnames(LoadedBackup{1}); % Get the list of field to process from one of the backup
for i = 1:length(fieldList)
    
    if isfield(LoadedBackup{2},fieldList(i))
        
        if ismember(fieldList(i), gridInfoListExt) % if grid information, simply copy of the field
            disp(['grid info : ' fieldList{i}]);

            if strcmp(fieldList{i},'Macrocaetes') % mean(A.macro, B.macro)
                disp(['mean : ' fieldList{i}]);
                DOA_backup.Macrocaetes = (LoadedBackup{1}.Macrocaetes + LoadedBackup{2}.Macrocaetes) ./ 2;    
            elseif strcmp(fieldList{i},'AreaRatios') % min(A.ar, B.ar)
                disp(['min : ' fieldList{i}]);
                eval(['DOA_backup.' fieldList{i} ' = min(LoadedBackup{1}.' fieldList{i} ',LoadedBackup{2}.' fieldList{i} ');']);
            elseif strncmp(fieldList{i},'error',5) % max(A.error, B.error)
                disp(['max : ' fieldList{i}]);
                eval(['DOA_backup.' fieldList{i} ' = max(LoadedBackup{1}.' fieldList{i} ',LoadedBackup{2}.' fieldList{i} ');']);
            else
                eval(['DOA_backup.' fieldList{i} ' = LoadedBackup{1}.' fieldList{i} ';']); % copied without change (1.5)
            end
            
        else
            if ~isempty(strfind(fieldList{i},'_Smap')) % if Significance map, logical AND(A.map, B.map)
                disp(['logical AND : ' fieldList{i}]);
                eval( ['DOA_backup.' fieldList{i} ' = or( LoadedBackup{1}.' fieldList{i} ', LoadedBackup{2}.' fieldList{i} ');'] );
            else
                disp(['Substract : ' fieldList{i}]);
                eval( ['DOA_backup.' fieldList{i} ' = LoadedBackup{1}.' fieldList{i} ' - LoadedBackup{2}.' fieldList{i} ';'] );
            end
        end
    end 
end

%%% save results
disp('Saving DOA backup ...')
if ~exist(OutputFolderName,'dir')
    mkdir(OutputFolderName);
end
save([OutputFolderName filesep 'DOA_backup'],'-struct', 'DOA_backup');


%% Plot %%

if ~isempty(Qname) && PLOT.plot
    
    %%% reload the backup (not necessary, we still do it as a verification process)
    GRID = load([OutputFolderName filesep 'DOA_backup']); % 1.5
    
    %%% check dimension (compatibility for single time point plot)
    Lt = size( eval(['GRID.' Qname{1}]) ,4);
    if Lt == 0
        Lt = size( eval(['GRID.' Qname{1}]) ,3);
    end
    
    % Extracting FULL or REG data for plot (1.6)
    TYPE = GRID.REG;
    if PLOT.makeItFull
        TYPE = GRID.FULL;
    end
    GRID.xywh = TYPE.xywh;
    GRID.centroids = TYPE.centroids;
    GRID.sizeImageX = TYPE.sizeImageX;
    GRID.sizeImageY = TYPE.sizeImageY;
    
    DISPLAY.yMid = TYPE.yMid;
    DISPLAY.macrocaetes = TYPE.macrocaetes;
    
    Qimage = ones(GRID.sizeImageY, GRID.sizeImageX);

    for t = 1:Lt
        
        % Sets the time and n variable (mod 1.4)
        timeStart = Time_str2dec(GRID.TimeArray{t,1}) + timeShift;                % turns time string into decimal and adds timeShift
        timeEnd = Time_str2dec(GRID.TimeArray{t,2}) + timeShift;
        timeStart = Time_dec2str(timeStart);                                      % switches back to string
        timeEnd = Time_dec2str(timeEnd);
        DISPLAY.time = [timeStart ' - ' timeEnd];
        DISPLAY.n = round( ( GRID.FrameArray(t,1) + GRID.FrameArray(t,2) ) / 2 );
        
        if Lt > 1 %(1.2)
            DISPLAY.step = t;
        end
        
        for j=1:length(Qname)
            
            if PLOT.Qparameters
                
                % Get the Qcolor, Qunits, Qsr, etc. corresponding to the Display (mod 1.5)
                [Pname,idx] = GetPname(Qname{j});
                if idx > 0
                    Qcolor  = eval(['allColors_' Pname '{idx}']);
                    Qunits  = eval(['allUnits_' Pname '{idx}']);
                    Qsr     = eval(['sr_' Pname '{idx}']);
                    Qsb     = eval(['srbar_' Pname '(idx)']);
                    QKillTr = eval(['killtrace_' Pname '(idx)']);
                end
            end
            
            % kill trace tag
            KillTrTag = '';
            if QKillTr
                KillTrTag = '_Tr=0';
            end
            
            % significant_map_wt = Significativity(D,dT_wt,t);
            if PLOT.significance
                temp_sign_map = eval( ['GRID.' Qname{j} '_Smap;'] );
                DISPLAY.significant_map = temp_sign_map(:,:,:,t);
            end
            
            % Plot image
            PlotField(Qname{j},QKillTr,Qcolor,Qunits,Qsr,Qsb,GRID,Qimage,DISPLAY);
            
            % Print plot
            if PLOT.print
                frameOutputFolder = [OutputFolderName filesep 'DOA_' QplotType fullTag]; % 1.5
                mkdir(frameOutputFolder)
                thisFilename = ['DOA_' DOAname '_' Qname{j} '_' timeStart 'to' timeEnd fullTag KillTrTag '_sr=' num2str(Qsr(1)) printFormat]; % 1.5
                if strcmp(imageExtension,'.svg')
                    plot2svg([frameOutputFolder filesep thisFilename],figure(1),'png');
                else
                    print (printFormat, printResolution, [frameOutputFolder filesep thisFilename]);
                end
            end
            close
        end % end for each quantities
    end % end for each time point
end


%% History %%

% 02/02/2017: 1.5 (Boris)
% - adjustments to make it work with AOA (2.7)

% 05/01/2017: 1.4 (Boris)
% - now uses "timeShift" to correct displayed time range of averaging
% - "Pname" became "PnameDOA"
% - "Qsrbar" became "Qsb"

% 19/03/2015: 1.3
% - code style modification and improvement
% - manage the significativity in boolean map

% 21/01/2015: 1.2
% - fixed bug when nb time point data == 1

% 13-01-2015 : 1.1
% - replot and print option
% - make difference between all the values in backups



