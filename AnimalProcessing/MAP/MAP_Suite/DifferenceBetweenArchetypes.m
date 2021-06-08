% DifferenceBetweenArchetypes (DBA)
% 
% Loads two AOA backups of two average animals and make the difference
% between the two.
%
% version 2.0 (DOA became DBA)
% Stephane Rigaud
% Boris Guirao

%% Creating DBA backup %%

BACKUP = MakeDBAbackup(LoadedBackupsRaw, DISPLAY, gridOverlap);

% Saves backup %%
disp('Saving DBA backup ...')
if ~exist(OutputPathName,'dir')
    mkdir(OutputPathName);
end
save([OutputPathName filesep 'DBA_backup'],'-struct', 'BACKUP');


%% Plot maps %%

if ~isempty(Qname) && PLOT.plot
    
    fprintf('Plotting DBA quantities ...\n')
    
    % Specific to DBA (the rest is copied from AOA)
    %---------------------------------------------------------
    Lt = size(BACKUP.AreaRatios,4);
    DISPLAY.Animal = DBAname;
%     BACKUP = load([OutputPathName filesep 'DBA_backup']); % added fullTag (2.6)
    %---------------------------------------------------------
    
    %%% Loading plot parameters (mod 2.9)
    % Extracting FULL or REG data for plot (2.7)
    TYPE = BACKUP.REG;
    if DISPLAY.makeItFull
        TYPE = BACKUP.FULL;
    end
    BACKUP.xywh = TYPE.xywh;
    BACKUP.Centroids = TYPE.Centroids;
    BACKUP.SizeImageX = TYPE.SizeImageX;
    BACKUP.SizeImageY = TYPE.SizeImageY;
    DISPLAY.yMid = TYPE.yMid;
    DISPLAY.macrocaetes = TYPE.Macrocaetes;
    
    Qimage = ones(BACKUP.SizeImageY, BACKUP.SizeImageX);
    
    %%% For each time point
    for t = 1:Lt
        
        % Correcting time display: ONLY relevant when TA quantities were added LAST in AOT backup(2.5)
        fixDtH = 0;
        timeShift = 0;
        % fixDtH = delta_t/60; % uncomment if 12h05 is displayed instead of 12h00 for instance
        
        % Set the time and n variable (mod 2.2)
        tmpTimeStart = TimeStr2Dec(BACKUP.TimeArray{t,1}) + timeShift - fixDtH;      % turns time string into decimal and adds timeShift, added fixDtH (2.5)
        tmpTimeEnd = TimeStr2Dec(BACKUP.TimeArray{t,2}) + timeShift;
        tmpTimeStart = TimeDec2Str(tmpTimeStart);                                            % switches back to string
        tmpTimeEnd = TimeDec2Str(tmpTimeEnd);
        DISPLAY.time = [tmpTimeStart ' - ' tmpTimeEnd];
        DISPLAY.n = round( ( BACKUP.FrameArray(t,1) + BACKUP.FrameArray(t,2) ) / 2 );
        
        if Lt > 1 % (1.2)
            DISPLAY.step = t;
        end
        
        %%% For each Display to be plot
        for j = 1:nQname
            
            if isfield(BACKUP,Qname{j}) % 2.3
                
                % Get the Qcolor, Qunits, Qsr, etc. corresponding to the Display
                [Pname,idx] = GetPname(Qname{j}); % now gets Pname (2.3)
                Qunits    = eval(['allUnits' Pname '{idx}']);
                Qcolor    = eval(['allColors' Pname '{idx}']);
                Qsr       = eval(['scaleRatio' Pname ]);
                Qscalebar = eval(['scaleBarLength' Pname ]);
                QKillTr   = eval(['killMeanTrace' Pname ]);
                Qsr = GetPlotScaleParameter(Qsr, idx);
                Qscalebar = GetPlotScaleParameter(Qscalebar, idx);
                QKillTr = GetPlotScaleParameter(QKillTr, idx);
                
                % Get the Qcolor, Qunits, Qsr, etc. corresponding to the Display
                [Pname,idx] = GetPname(Qname{j}); % now gets Pname (2.3)
                Qunits    = eval(['allUnits' Pname '{idx}']);
                Qcolor    = eval(['allColors' Pname '{idx}']);
                
                % tag indicating mean Tr = 0 for naming files
                KillTrTag = '';
                if QKillTr
                    KillTrTag = '_Tr=0';
                end
                
                % Load significance map
                if PLOT.significance && ~singleAnimal
                    DISPLAY.SignificanceMap = eval(['BACKUP.' Qname{j} '_Smap(:,:,:,t)']);
                end
                
                % Plot image
                PlotField(Qname{j}, QKillTr, Qcolor , Qunits , Qsr, Qscalebar, BACKUP, Qimage, DISPLAY);
                
                % Plot print (removed "Frames_" prefix and added fullTag in filenames)(2.6)
                thisFilename = ['DBA_' Qname{j} '_' tmpTimeStart 'to' tmpTimeEnd fullTag KillTrTag '_sr=' num2str(Qsr(1)) printFormat];                
%                 thisFilename = ['AOA_' Qname{j} '_' tmpTimeStart 'to' tmpTimeEnd KillTrTag fullTag '_sr=' num2str(Qsr(1)) imageExtension]; % 1.7, added fullTag (2.5), use timeStart/End (2.7)
                plotFolder = [OutputPathName filesep 'DBA_' plotType fullTag];
                if ~exist(plotFolder, 'dir')
                    mkdir(plotFolder);
                end
                if strcmp(imageExtension,'.svg')
                    plot2svg([plotFolder filesep thisFilename],figure(1),'png');
                else
                    print(printFormat, printResolution, [plotFolder filesep thisFilename]);
                end
                close
                
            else % 2.3
                disp(['DBA WARNING: quantity "' Qname{j} '" was not found in "BACKUP"!']);
                disp('Please make sure that it is the name assigned to this quantity AND that it was averaged over time and stored in each animal AOT backups.')
            end
        end % end for each quantities
    end % end for each time point
end % end if plot

%% History %%

% 03/06/2020: 2.0 (Boris) DOA became DBA
% - now using function "MakeDBAbackup" to create BACKUP

% 07/02/2017: 1.6 (Boris)
% - removed "DOAname" from file name of images (was too long)

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



