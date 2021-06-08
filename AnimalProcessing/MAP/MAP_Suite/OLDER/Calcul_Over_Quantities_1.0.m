% Loads one backup and make an operation between a quantity and a list of quantities
%
% Stephane Rigaud
% version 1
%

%% process

%for i = length(Qname1)
if strcmp(CalculType,'minus')
    disp(['Processing : ' QnameA ' - ' QnameB]);
    eval(['BACKUP.' QnameA CalculType QnameB ' = BACKUP.' QnameA ' - BACKUP.' QnameB ';']);
    if ~isempty(Oname)
        eval(['BACKUP.' QnameA CalculType QnameB '_Smap = or(BACKUP.' QnameA '_Smap, BACKUP.' QnameB '_Smap);']);
    end
end
if strcmp(CalculType,'plus')
    disp(['Processing : ' QnameA ' + ' QnameB]);
    eval(['BACKUP.' QnameA CalculType QnameB ' = BACKUP.' QnameA ' + BACKUP.' QnameB ';']);
    if ~isempty(Oname)
        eval(['BACKUP.' QnameA CalculType QnameB '_Smap = or(BACKUP.' QnameA '_Smap, BACKUP.' QnameB '_Smap);']);
    end
end
%end


%%% Save AOA mean structure
disp('Appending structure ...')
Previous = load([OutputPath filesep backupPathName]);
eval([backupName ' = catstruct(Previous,BACKUP);']);
save([OutputPath filesep backupPathName], '-struct', backupName, '-append');


%% plot
if PLOT.plot
    
    disp('Starting Plot ...')
    
    %%% initialisation of plot
    BACKUP = load([OutputPath filesep backupPathName]);
    
    Lt = size(BACKUP.FrameArray,1);
    
    [Pname,idx] = GetPname(QnameA);
    Qsr       = eval(['sr_' Pname '{idx}']);
    Qscalebar = eval(['srbar_' Pname '(idx)']);
    Qunits    = eval(['allUnits_' Pname '{idx}']);
    Qcolor    = eval(['allColors_' Pname '{idx}']);
    QKillTr   = eval(['killtrace_' Pname '(idx)']);
    
    KillTrTag = '';
    if QKillTr
        KillTrTag = '_Tr=0';
    end
    
    %%% get the display info parameters
    DISPLAY.Animal             = Animal;
    DISPLAY.plotType           = QplotType;
    DISPLAY.minimalInfoDisplay = minimalInfoDisplay;
    DISPLAY.minAEV             = minAEV;
    DISPLAY.scaleBarWidth      = scaleBarWidth;
    DISPLAY.gridDisplay        = gridDisplay;
    DISPLAY.lineWidth          = lineWidth;
    DISPLAY.pointSize          = pointSize;
    DISPLAY.EVstyles           = EVstyles;
    DISPLAY.signOpacities      = signOpacities;
    DISPLAY.fontSizeInfo       = fontSize;
    DISPLAY.imageFading        = imageFading;
    DISPLAY.errorPsMin         = 1;
    DISPLAY.errorDnPsMin       = 1;
    DISPLAY.errorFontSize      = 5;
    if PLOT.macrocaetes
        DISPLAY.macrocaetes = BACKUP.Macrocaetes;
    end
    
    
    %%% BG image
    size_image_x = PLOT.boxSize(1) * (BACKUP.size(2) * BACKUP.overlap + 1);
    size_image_y = PLOT.boxSize(2) * (BACKUP.size(1) * BACKUP.overlap + 1);
    if ~BACKUP.overlap
        size_image_x = PLOT.boxSize(1) * BACKUP.size(2);
        size_image_y = PLOT.boxSize(2) * BACKUP.size(1);
    end
    Qimage = ones([size_image_y , size_image_x]);
    
    
    %%% plot process
    for t = 1:Lt
        DISPLAY.time = [ BACKUP.TimeArray{t,1} ' - ' BACKUP.TimeArray{t,2} ];
        DISPLAY.n = round( ( BACKUP.FrameArray(t,1) + BACKUP.FrameArray(t,2) ) / 2 );
        
        if Lt > 1
            DISPLAY.step = t;
        end
        
        PlotQname = [QnameA CalculType QnameB];
        if ~isempty(Oname)
            eval(['BACKUP.AreaRatios = BACKUP.AreaRatios_' Pname ';']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(Oname)
            framePlot = BACKUP.FrameArray(t,2);
            pathbackup = [RAW.gridAnimalFolder filesep 'Backups' filesep Pname '_LGrid_' Animal '_' num2str(framePlot,['%0' num2str(RAW.nDigits) 'd']) '.mat'];
            BU = load(pathbackup);
            iGRID = BU.(['GRID_' Pname]);
            DISPLAY.Lcentroids = iGRID.Lcentroids;
            DISPLAY.contour_indices = iGRID.contour_indices;
            Qimage = imread([RAW.rawPathFolder filesep 'Output_results' filesep RAW.filename num2str(framePlot, ['%0' num2str(RAW.nDigits) 'd']) '.' RAW.imageFormat]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Significance
        if PLOT.significance
            temp_sign_map = eval(['BACKUP.' PlotQname '_Smap']);
            DISPLAY.significant_map = temp_sign_map(:,:,:,t);
        end
        
        PlotField(PlotQname,QKillTr,Qcolor,Qunits,Qsr,Qscalebar,BACKUP,Qimage,DISPLAY);
        
        %%% print
        if PLOT.print
            this_filename = [Oname '_' PlotQname '_' BACKUP.TimeArray{t,1} 'to' BACKUP.TimeArray{t,2} KillTrTag '_sr=' num2str(Qsr) imageExtension]; % 1.7
            mkdir([OutputPath filesep 'Frames_COQ_' QplotType]);
            if strcmp(imageExtension,'.svg')
                plot2svg([OutputPath filesep 'Frames_COQ_' QplotType filesep this_filename],figure(1),'png');
            else
                print(printFormat, printResolution, [OutputPath filesep 'Frames_COQ_' QplotType filesep this_filename]);
            end
        end
        close
        
    end
end

