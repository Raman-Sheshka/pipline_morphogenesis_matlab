

Version = 4.0;
% Stephane Rigaud
% Boris Guirao






%% Backup initialisation %%

newBackup = struct();
gridInfoListExt = [gridInfoList ; 'AreaRatios'];
for i=1:size(gridInfoListExt,1)
    if isfield(aBackup,gridInfoListExt{i})
        eval(['newBackup.' gridInfoListExt{i}  '= aBackup.' gridInfoListExt{i} ';']);
    end
end




%% Initialisation of the unitary tensors %%

[by, bx, bz, bt, ba] = size(Qproj);

thetaB = GetAngleMap(Qproj); % return the angle map in rad
sigma0 = repmat(reshape([1 0 ; 0  1], [1 1 4 1 1]), [by bx 1 bt ba]);
sigma1 = repmat(reshape([0 1 ; 1  0], [1 1 4 1 1]), [by bx 1 bt ba]);
sigma3 = repmat(reshape([1 0 ; 0 -1], [1 1 4 1 1]), [by bx 1 bt ba]);

sigmaIso  = sigma0;
sigmaDev  = repmat( cos(thetaB.*2), [1 1 4 1 1]) .* sigma3 + repmat(sin(thetaB.*2), [1 1 4 1 1]) .* sigma1;
sigmaDevO = repmat(-sin(thetaB.*2), [1 1 4 1 1]) .* sigma3 + repmat(cos(thetaB.*2), [1 1 4 1 1]) .* sigma1;

unitaryTensorList = {sigmaIso ; sigmaDev ; sigmaDevO};
tagUnitaryTensor  = {   'i'   ;    'd'   ;   'do'   };




%% Projection of quantities over uProjTensor %%

for i = 1:nQname
    for t = 1:length(tagUnitaryTensor)
        
        fprintf(['\tComputing ' Qname{i} '.' tagUnitaryTensor{t} '(' uQname ') ...']);
        
        if isfield(aBackup,Qname{i})
            
            currentQ = eval(['aBackup.' Qname{i} ';']);
            uProjTensor = unitaryTensorList{t};
            
            if size(currentQ, 3) == 4
                
                %%% managing quantities size
                if size(currentQ, 4) > size(uProjTensor, 4)
                    uProjTensor = repmat(uProjTensor, [1 1 1 size(currentQ, 4) 1]);
                elseif size(currentQ, 4) < size(uProjTensor, 4)
                    currentQ = repmat(currentQ, [1 1 1 size(uProjTensor, 4) 1]);
                end
                
                %%% compute projection
                AdotB = ScalarProductFunction(currentQ, uProjTensor);
                
                %%% manage significance map
                eval(['newBackup.' Qname{i} 'dot' uQname '_' tagUnitaryTensor{t} ' = AdotB;']);
                if ~singleAnimal % 3.7
                    if strcmp(tagUnitaryTensor{t},'d') || strcmp(tagUnitaryTensor{t},'do')
                        
                        temp_map = eval(['aBackup.' uQname '_Smap(:,:,2,:,:)']);
                        B_Smap = cat( 3, temp_map, temp_map );
                        eval(['newBackup.' Qname{i} 'dot' uQname '_' tagUnitaryTensor{t} '_Smap = B_Smap;']);
                    elseif strcmp(tagUnitaryTensor{t},'i')
                        
                        temp_map = eval(['aBackup.' uQname '_Smap(:,:,1,:,:)']);
                        B_Smap = cat( 3, temp_map, temp_map );
                        eval(['newBackup.' Qname{i} 'dot' uQname '_' tagUnitaryTensor{t} '_Smap = B_Smap;']);
                    end
                end
                fprintf(' Done\n');
            else
                fprintf(' Skipped\n');
            end
        else
            fprintf(' Skipped\n');
        end
    end
end

%%% save, or update if already existing, the new backup
if ~exist(OutputPathName,'dir')
    mkdir(OutputPathName);
end

thisFilename = [OutputPathName filesep 'POA_' uAnimal '_' tagProjectionTime '.mat'];
if ~exist(thisFilename,'file')
    disp('Saving new structure ...')
    save(thisFilename,'-struct', 'newBackup');
else
    disp('Appending structure ...')
    save(thisFilename,'-struct', 'newBackup', '-append');
end




%% Plot %%

if ~isempty(Qname) && PLOT.plot && ~skipPOAplots
    fprintf('Plotting projected quantities ...\n')
    
    %%% Load backup
    BACKUP = load([OutputPathName filesep 'POA_' uAnimal '_' tagProjectionTime]); % added fullTag (2.6)
    
    Lt = size(BACKUP.FrameArray,1);
    
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
    
    
    %%% plot process
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
            
            for p = 1:length(tagProjPlot)
                
                PlotQname = [Qname{j} 'dot' uQname '_' tagProjPlot{p}];
                
                if isfield(BACKUP,PlotQname)
                    
                    POAname = ['POA_' Qname{j} '.' uQname '_' tagProjPlot{p} '_'  uAnimal  '_' tagProjectionTime];
                    DISPLAY.Animal = POAname;
                    [Pname,idx] = GetPname(Qname{j});
                    Qcolor  = eval(['allColors' Pname '{idx}']);
                    QKillTr = eval(['killMeanTrace' Pname ]);
                    if length(QKillTr)>1
                        QKillTr = QKillTr(idx);
                    end
                    
                    Qunits = eval(['allUnits' Pname '{idx}']);
                    Qsr    = eval(['scaleRatio' Pname '{idx}']) .* scaleRatioFactor;
                    Qsb    = eval(['scaleBarLength' Pname '(idx)']);
                    
                    %%% tag indicating mean Tr = 0 for naming files
                    KillTrTag = '';
                    if QKillTr
                        KillTrTag = '_Tr=0';
                    end
                    
                    La = size(eval(['BACKUP.' PlotQname]),5);
                    if La > 1 % (1.2)
                        DISPLAY.animalIdx = t;
                    end
                    
                    %%% Significance
                    if PLOT.significance
                        DISPLAY.significant_map = eval(['BACKUP.' PlotQname '_Smap(:,:,:,t)']);
                    end
                    
                    %%% Plot value as a surface instead of radius
                    S = sign( eval(['BACKUP.' PlotQname ';']) );
                    V = sqrt( abs( eval(['BACKUP.' PlotQname ';']) ) );
                    eval(['BACKUP.' PlotQname ' = S .* V;']);
                    
                    PlotField(PlotQname, QKillTr, Qcolor, Qunits, Qsr, Qsb, BACKUP, Qimage, DISPLAY);
                    
                    % Plot print (removed "Frames_" prefix and added fullTag in filenames)(2.6)
                    thisFilename = ['POA_' PlotQname '_' tmpTimeStart 'to' tmpTimeEnd KillTrTag fullTag '_sr=' num2str(Qsr(1)) imageExtension]; % 1.7, added fullTag (2.5), use timeStart/End (2.7)
                    plotFolder = [OutputPathName filesep 'POA_' plotType fullTag];
                    if ~exist(plotFolder, 'dir')
                        mkdir(plotFolder);
                    end
                    if strcmp(imageExtension,'.svg')
                        plot2svg([plotFolder filesep thisFilename], figure(1), 'png');
                    else
                        print(printFormat, printResolution, [plotFolder filesep thisFilename]);
                    end
                    close
                end
            end
        end 
    end
end


%% History %%

% 13/06/2018: 4.0 (Stephane)
% - update for new SAP and MAP script
% - add a fifth dimension to compute projection on average animal (A1) and per animal (A2 to An)
% - remove single animal processing (single animal only done in SAP now)

% 06/02/2017: 3.7 (Boris)
% - adjustments to process single animals
% - now directly compares first 2 characters of A to the list of scalar quantities to skip A projection
% - STOPPED saving non-tensor quantities for LTA plots: now LTA loads AOA or DOA AND POA backups simultaneously
% - changed filenames of POA images being saved

% 02/02/2017: 3.6 (Boris)
% - look for "3.6" to see changes

% 27/01/2017: 3.5 (Boris)
% - now still runs when quantities listed in Qname are not available (or are non-tensors cf 3.4)
% - now saving non-tensor quantities (such as dnA) and their Smap UNTOUCHED in POA backups for convenience in LTA plot

% 26/01/2017: 3.4 (Boris)
% - now runs even when non-tensor quantites are listed in Qname in AIA_MultiOperation (they're just skipped)
% - now takes into account "timeShift" when naming files being saved

% 23/01/2017: 3.3 (Boris)
% - adjustments for compatibility with AOA v2.7 and the use of unique backups for regular or full case (makeItFull = true)
% - accordingly loads structures REG and FULL stored in AOA backups

% 20/01/2017: 3.2 (Boris)
% - ONLY using AOT backups now!!
% - moved most filling of DISPLAY to AOA_MultiOperation

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