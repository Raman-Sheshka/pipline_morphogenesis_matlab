% ProjectionOverAnimals (POA)
%
% Is called either by:
% AIA_Info_Multi to run on AOA or DOA data
% X_AIA_Info to run on the animal X
%
% LoadedBackup contain the input backup structure containing the quantitites to be contracted
% Qname is a list of two string which are the name of the quantitites to be contracted
% deviator a flag to say if we process only the deviatoric part
% denominateur an option to devide the contradted product by a value (currently only the squared norm of B)
%
% print ; plot ; replot are flag to manage different part of the process
%
% version 3.1
% Stephane Rigaud
% Boris Guirao

%% process


%%% Backup initialisation
newBackup = struct();
GridInfo = {'xywh' ; 'size' ; 'Lcentroids' ; 'overlap' ; 'color' ; 'centroids' ; 'ULCs' ; 'lineWidth' ; 'TimeArray' ; 'FrameArray' ; 'RConds' ; 'Macrocaetes' ; 'AreaRatios_TA' ; 'AreaRatios_SM' ; 'AreaRatios_AOS' ; 'AreaRatios_VMM'};
for i=1:size(GridInfo,1)
    if isfield(BACKUP,GridInfo{i})
        eval(['newBackup.' GridInfo{i}  '= BACKUP.' GridInfo{i} ';']);
    end
end
% single individual case
% if isempty(OnamePOA)
%     eval(['[wy,wx,wt] = size(BACKUP.AreaRatios_' Pname ');']);
%     if strcmp(Pname,'TA')
%         reshapeR = reshape(BACKUP.RConds,[wy,wx,1,wt]);
%         newBackup.RConds = reshapeR;
%     end
%     eval(['reshapeW = reshape(BACKUP.AreaRatios_' Pname ',[wy,wx,1,wt]);']);
%     eval(['newBackup.AreaRatios_' Pname ' = reshapeW;']);
% end

%%% Initialisation of the projection unitary tensor
[by, bx, bz, bt] = size(B);

thetaB = GetAngleMap(B); % return the angle map in rad
sigma0 = repmat(reshape([1 0 ; 0  1],[1 1 4 1]),[by bx 1 bt]); 
sigma1 = repmat(reshape([0 1 ; 1  0],[1 1 4 1]),[by bx 1 bt]); 
sigma3 = repmat(reshape([1 0 ; 0 -1],[1 1 4 1]),[by bx 1 bt]); 

sigmaIso  = sigma0;
sigmaDev  = repmat( cos(thetaB.*2),[1 1 4 1]) .* sigma3 + repmat(sin(thetaB.*2),[1 1 4 1]) .* sigma1;
sigmaDevO = repmat(-sin(thetaB.*2),[1 1 4 1]) .* sigma3 + repmat(cos(thetaB.*2),[1 1 4 1]) .* sigma1;

unitaryMatrices = {sigmaIso ; sigmaDev ; sigmaDevO};

%%% Projection of A over uB
for t = 1:length(tagTensor)
    
    disp(['- tagTensor ' tagTensor{t}]);
    Bu = unitaryMatrices{t};
        
    for i = 1:length(Qname)
        disp(['    Processing ' Qname{i} '.' tagTensor{t} '(' uQname ') ...']);
        
        % loading A
        A = eval(['BACKUP.' Qname{i} ';']);
        
        sizedBu = Bu;
        if size(A,4) ~= size(Bu,4)
            sizedBu = repmat(Bu,[1 1 1 size(A,4)]);
        end
        
        % projection function
        POA = ScalarProductFunction(A,sizedBu);
        
        % manage significance map
        if ~isempty(OnamePOA)
            if strcmp(tagTensor{t},'d') || strcmp(tagTensor{t},'do')
                temp_map = eval(['BACKUP.' uQname '_Smap(:,:,2,:)']);
                B_Smap = cat( 3, temp_map, temp_map );
                eval(['newBackup.' Qname{i} 'dot' uQname '_' tagTensor{t} '_Smap = B_Smap;']);
            elseif strcmp(tagTensor{t},'i')
                temp_map = eval(['BACKUP.' uQname '_Smap(:,:,1,:)']);
                B_Smap = cat( 3, temp_map, temp_map );
                eval(['newBackup.' Qname{i} 'dot' uQname '_' tagTensor{t} '_Smap = B_Smap;']);
            end
        end
        eval(['newBackup.' Qname{i} 'dot' uQname '_' tagTensor{t} ' = POA;']);
        
    end
end

% save backup under template POA_uAnimal_uTimeWidth_uTimeOverlap
% the projected animal is already specified in the folder name
eval(['POA_' uAnimal '_' tagProjectionTime ' = newBackup;']);

% save, or update if already existing, the new backup
if ~exist(outputPath,'dir'); mkdir(outputPath); end
if ~exist([outputPath filesep 'POA_' uAnimal '_' tagProjectionTime '.mat'],'file')
    disp('Saving new structure ...')
    save([outputPath filesep 'POA_' uAnimal '_' tagProjectionTime ],'-struct', ['POA_' uAnimal '_' tagProjectionTime]);
else
    disp('Appending structure ...')
    Previous = load([outputPath filesep 'POA_' uAnimal '_' tagProjectionTime]);
    eval(['POA_' uAnimal '_' tagProjectionTime ' = catstruct(Previous,POA_' uAnimal '_' tagProjectionTime ');']);
    save([outputPath filesep 'POA_' uAnimal '_' tagProjectionTime],'-struct', ['POA_' uAnimal '_' tagProjectionTime], '-append');
end


%% plot
if PLOT.plot
    
    disp('Starting Plot ...')
    
    %%% initialisation of plot
    BACKUP = load([outputPath filesep 'POA_' uAnimal '_' tagProjectionTime ]);
    Lt = size(BACKUP.FrameArray,1);
    
    %%% plot parameters
    DISPLAY.plotType = QplotType;
    DISPLAY.minimalInfoDisplay = minimalInfoDisplay;
%     DISPLAY.minAEV = minAEV;
    DISPLAY.scaleBarWidth = scaleBarWidth;
    DISPLAY.gridDisplay = gridDisplay;
    DISPLAY.lineWidth = lineWidth;
    DISPLAY.pointSize = pointSize;
    DISPLAY.EVstyles = EVstyles;
    DISPLAY.signOpacities = signOpacities;
    DISPLAY.fontSizeInfo = fontSize;
    DISPLAY.imageFading = imageFading;
    DISPLAY.errorPsMin = 1;
    DISPLAY.errorDnPsMin = 1;
    DISPLAY.errorFontSize = 5;
    if isfield(BACKUP,'Macrocaetes');
        DISPLAY.macrocaetes = BACKUP.Macrocaetes;
    end
    
    %%% BG image
    size_image_x = PLOT.boxSize(1) * (BACKUP.size(2) * BACKUP.overlap + 1);
    size_image_y = PLOT.boxSize(2) * (BACKUP.size(1) * BACKUP.overlap + 1);
    if BACKUP.overlap == 0
        size_image_x = PLOT.boxSize(1) * BACKUP.size(2);
        size_image_y = PLOT.boxSize(2) * BACKUP.size(1);
    end
    Qimage = ones([size_image_y , size_image_x]);
    % Qimage = ones([2016 , 3559]);
    
    
    
    %%% plot process
    for t=1:Lt
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(OnamePOA)
            framePlot = BACKUP.FrameArray(t,2);
            pathbackup = [RAW.gridAnimalFolder filesep 'Backups' filesep Pname '_LGrid_' Animal '_' num2str(framePlot,['%0' num2str(RAW.nDigits) 'd']) '.mat'];
            BU = load(pathbackup);
            iGRID = BU.(['GRID_' Pname]);
            DISPLAY.Lcentroids = iGRID.Lcentroids;
            DISPLAY.contour_indices = iGRID.contour_indices;
            Qimage = imread([RAW.rawPathFolder filesep 'Output_results' filesep RAW.filename num2str(framePlot, ['%0' num2str(RAW.nDigits) 'd']) '.' RAW.imageFormat]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Sets the time and n variable (mod 3.1)
        timeStart = Time_str2dec(BACKUP.TimeArray{t,1}) + timeShift;                % turns time string into decimal and adds timeShift
        timeEnd = Time_str2dec(BACKUP.TimeArray{t,2}) + timeShift;
        timeStart = Time_dec2str(timeStart);                                      % switches back to string
        timeEnd = Time_dec2str(timeEnd);
        DISPLAY.time = [timeStart ' - ' timeEnd];
%         DISPLAY.time = [ BACKUP.TimeArray{t,1} ' - ' BACKUP.TimeArray{t,2} ];
        DISPLAY.n = round( ( BACKUP.FrameArray(t,1) + BACKUP.FrameArray(t,2) ) / 2 );
        
        if Lt > 1
            DISPLAY.step = t;
        end
        
        for n = 1:length(Qname)
            
            PlotQname = [Qname{n} 'dot' uQname '_' tagPlot];
            POAname = ['POA_' uAnimal '_' Qname{n} '.' uQname '_' tagPlot ];
            DISPLAY.Animal = POAname;
            
            eval(['BACKUP.AreaRatios = min( BACKUP.AreaRatios_' Pname ', BACKUP.AreaRatios_' uPname ');']);

            [~,idx] = GetPname(Qname{n});
            Qcolor  = eval(['allColors_' Pname '{idx}']);
            %QunitsPOA  = eval(['allUnits_' Pname '{idx}']);
            %QsrPOA     = eval(['sr_' Pname '{idx}']);
            %QsbPOA  = eval(['srbar_' Pname '(idx)']);
            QKillTr = eval(['killtrace_' Pname '(idx)']);
            
            %%% Significance
            if PLOT.significance
                temp_sign_map = eval(['BACKUP.' PlotQname '_Smap']);
                DISPLAY.significant_map = temp_sign_map(:,:,:,t);
            end
            
            %%% Plot value as a surface instead of radius
            S = sign( eval(['BACKUP.' PlotQname ';']) );
            V = sqrt( abs( eval(['BACKUP.' PlotQname ';']) ) );
            eval(['BACKUP.' PlotQname ' = S .* V;']);
            
            PlotField(PlotQname,QKillTr,Qcolor,QunitsPOA,QsrPOA,QsbPOA,BACKUP,Qimage,DISPLAY);
            
            %%% print
            if PLOT.print
                this_filename = ['POA_' uAnimal '_' PlotQname '_' tagProjectionTime '_' BACKUP.TimeArray{t,1} 'to' BACKUP.TimeArray{t,2} '_sr=' num2str(Qsr) imageExtension]; % 1.7
                mkdir([outputPath filesep 'Frames_POA_' QplotType]);
                if strcmp(imageExtension,'.svg')
                    plot2svg([outputPath filesep 'Frames_POA_' QplotType filesep this_filename],figure(1),'png');
                else
                    print(printFormat, printResolution, [outputPath filesep 'Frames_POA_' QplotType filesep this_filename]);
                end
            end
            close
            
        end
    end
end


%% Historique

% 04/01/2017: 3.1 (Boris)
% - use of "OnamePOA" instead of "Oname"
% - now uses "timeShift" to correct displayed time range of averaging
% - Qunits, Qsr, Qsb became QunitsPOA, QsrPOA, QsbPOA

% 07/08/2015: v3
% - projection on unitary tensor sigma0 sigmaDev sigmaOrtho
% - use of the new TensorData

% 21/07/2015: v2.5
% - remove normalisation
% - add global or instant contraction

% 19/03/2015: v2.0
% - Major modification of the code
