% LocalTimeAnalysis
% Merge all MAP backups and allow the user to plot, in a specific region, one or more quantities over time

version = 4.0;
% Stephane Rigaud
% Boris Guirao

%% Code %%

%%% Backup size check with first quantity listed in "QnameLTA" (revamped 4.0)
Quantity1Name = [QnameLTA{1} Qtags{1}];
Quantity1NameAlt = QnameLTA{1};         % without "iso"
if isfield(BACKUP, Quantity1Name) 
    [Ly, Lx, ~, Lt, La] = size(BACKUP.(Quantity1Name)); 
else
    [Ly, Lx, ~, Lt, La] = size(BACKUP.(Quantity1NameAlt));
end

TimeDecArray = mean(cellfun(@TimeStr2Dec,BACKUP.TimeArray),2);
startAPF = round(TimeDecArray(1));
endAPF   = round(TimeDecArray(end));

%%% Define Areas to process
% shift grid based on origin
gridOlap = 1 ./ (1-gridOverlap);
oriBox = cellfun(@(x)sum(abs(x - BACKUP.REG.xywh(1:2))), BACKUP.ULCs, 'UniformOutput', false);
[bOriY,bOriX] = find(cell2mat(oriBox) == 0);
ShiftBoxCoord = cellfun(@(x)(floor(x .* gridOlap) - [bOriX bOriY]),BACKUP.Coordinates,'UniformOutput',false);

%NB: Not the case currently but end goal is that, whatever the mean, we get
% cell array of origin based coordinate of boxes.
[res] = InputParserLTA(Qboxes);
switch res
    case 'gui'
        % graphic interface for area selection
        [~,Qboxes] = GuiLTA([bOriX bOriY], BACKUP, clickTime, PLOT);
        boxTag = ['(' regexprep(num2str(Qboxes{1}), ' *' , ',' ) ') ('  regexprep(num2str(Qboxes{end}), ' *' , ',' ) ')'];
        
    case 'key'
        % toBeExtended by other keyword
        Qboxes = ShiftBoxCoord(:);
        boxTag = ['(full)'];
        
    case 'mb'
        boxTag = ['(' regexprep(num2str(Qboxes{1}), ' *' , ',' ) ') ('  regexprep(num2str(Qboxes{end}), ' *' , ',' ) ')'];
        
    case 'ma'
        point1 = Qboxes{1}(1:2) + [bOriX bOriY];
        point2 = Qboxes{1}(3:4) + [bOriX bOriY];
        [listIdx] = ComputeAreaCoord(point1, point2, [Ly Lx]);
        [yl, xl] = ind2sub(GUIgrid.Size, listIdx);
        Qboxes = mat2cell([xl yl], ones(numel(xl),1), 2);
        Qboxes = cellfun(@(x)(x - ([bOriX bOriY])), Qboxes, 'UniformOutput', false);
        boxTag = ['(' regexprep(num2str(Qboxes{1}), ' *' , ',' ) ') ('  regexprep(num2str(Qboxes{end}), ' *' , ',' ) ')'];
        
    otherwise
        return
end

%%% shift back boxes coordinate and switch them into indexes
QboxesShift = cellfun(@(x)(x + ([bOriX bOriY])), Qboxes, 'UniformOutput', false); % reassigns [1 1] to top left grid compartment
PointList = cell2mat(QboxesShift); % removed the +1 that was a mistake (4.0)
% PointList = cell2mat(Qboxes) + 1;
idxList = sub2ind([Ly Lx], PointList(:,2), PointList(:,1));

%%% reshape boxes into logical map for data extraction
nbBoxes = numel(Qboxes);
BoxIndexMap = false(Ly, Lx);
BoxIndexMap(idxList) = true;
reshapeBoxIndexMap = repmat(BoxIndexMap, [1 1 1 Lt La]);
boxTagFile = regexprep(boxTag,' *', '');

%%% iteration per tags
for t = 1:length(Qtags)
    
    Total = zeros(1,Lt,La);
    Qmean = cell(1,size(QnameLTA,1));
    legendQname = {}; % 3.5
    
    %%% iteration per quantities
    for q = 1:length(QnameLTA)
        
        % Get information the data
        QuantityName = [QnameLTA{q} Qtags{t}];
        QuantityNameAlt = QnameLTA{q};
        if t == 1 && ~isfield(BACKUP, QuantityName) && isfield(BACKUP, QuantityNameAlt) % tag "iso" or nothing (4.0)
            QuantityName = QuantityNameAlt; % for scalar quantities other than iso
        end
        [Q1,Q2] = QuantityNameSpliter(QuantityName);
        Pname = GetPname(Q1);
        uPname = '';
        if ~isempty(Q2)
            uPname = GetPname(Q2);
        end
        
        % weight calculation (mod 3.4)
        subAreaRatios = reshape(BACKUP.AreaRatios(reshapeBoxIndexMap),[1 nbBoxes 1 Lt La]);
        Weight = subAreaRatios .^ 2;
        
        % Calculate qunatity to plot
        if isfield(BACKUP, QuantityName) % 3.5
            fprintf(['\t\tExtracting ' QnameLTA{q} Qtags{t} ' ...'])
            
            % extract quantity of the crop zone
            Quantity = eval(['BACKUP.' QuantityName ';']);
            if size(Quantity,3) ~= 1 % not a scalar
                Quantity = TensorNormMap(Quantity);
            end
            subQuantity = reshape(Quantity(reshapeBoxIndexMap),[1 nbBoxes 1 Lt La]);
            
            % Calculate weighted mean of the crop zone
            weightedSubQuantity = subQuantity .* Weight;
            weightedSubMean = nansum(weightedSubQuantity,2) ./ sum(Weight,2);
            weightedSubMean = reshape(weightedSubMean,[1 Lt La]);
            
            % save
            Qmean{q} = weightedSubMean;
            Total = Total + weightedSubMean;
            %--------------------------------------------------------------
            fprintf(' Done\n')
        end
    end
    
    
    %% Plot  %%
    
    % Updating "QmaxValues" (4.0)
    if ~isempty(QmaxValues)
        if length(QmaxValues) == length(QnameLTA) % right number of maxValues was specified
            QmaxValuesMod = QmaxValues;
        elseif length(QmaxValues) == 1
            QmaxValuesMod = repmat(QmaxValues{1},length(QnameLTA),1); % repeating unique value specified
            QmaxValuesMod = num2cell(QmaxValuesMod);
        else
            warndlg('ERROR: parameter "QmaxValues" must either contain, none, 1, or the same number of entries as "QnameLTA"!','LTA ERROR!!')
            return
        end
    end
    
    sth2PlotTF = any(~cellfun(@isempty, Qmean)); % 4.0

    if sth2PlotTF % only plots if there is something to plot (4.0)

        fprintf('\t\tPlotting values ...')
        
        % define figure for plot
        figure(1);
        box on % 4.0
        if ~isempty(Qrange)
            if length(Qrange) == 1
                Qrange = repmat(Qrange,4,1); % repeats range 4 times to do '', 'd', 'i', 'do' plots (3.5)
            end
            axis([ startAPF endAPF Qrange{t}]);
        end
        hold on
        
        
        for q = 1:length(QnameLTA)
            
            QuantityName = [QnameLTA{q} Qtags{t}];
            QuantityNameAlt = QnameLTA{q};
            if t == 1 && ~isfield(BACKUP, QuantityName) && isfield(BACKUP, QuantityNameAlt) % tag "iso" or nothing (4.0)
                QuantityName = QuantityNameAlt; % for scalar quantities other than "iso" part
            end
            
            if isfield(BACKUP, QuantityName) % 3.5
                
                % manage complex quantity name
                [Q1,Q2] = QuantityNameSpliter(QuantityName);
                plotName = QuantityName;
                [Pname,idx] = GetPname(Q1);
                uPname = '';
                if ~isempty(Q2)
                    uPname = GetPname(Q2);
                end
                Qcolor = eval(['allColors' Pname '{idx}']);
                Qunits = eval(['allUnits'  Pname '{idx}']);
                %--------------------------------------------------------------
                % get value and mean value
                value = Qmean{q};
                value = Qfactor(q)*value; % 3.4
                if length(value) == 1
                    Lt = 2;
                    value = Qmean{q} .* ones( size( [startAPF:(endAPF-startAPF)/(Lt-1):endAPF] ) );
                end
                % mod 4.0
                if strcmp(QplotType,'cum')
                    value = nancumsum(value*dtH);
                    if strcmp(Qunits,'h^{-1}')
                        Qunits = '';
                    else
                        Qunits = [Qunits '.h'];
                    end
                end
                %--------------------------------------------------------------
                
                if length(QnameLTA) == 1
                    %%% Plotting single value but with multianimal details
                    legendAnimal = cat(1,'mean Animal',avgAnimals);
                    CM = spring(numel(legendAnimal) -1); % +1 is for the reference movie
                    for a = 1:length(legendAnimal)
                        %--------------------------------------------------------------
                        % plot, legend, color, units
                        if a == 1
                            QLineWidth = 2;
                            QLineStyle = '-';
                        else
                            Qcolor = rand(1,3);
                            QLineWidth = 1;
                            QLineStyle = '--';
                        end
                        
                        % 4.0
                        if strcmp(QplotRenorm,'renorm')
                            if isempty(QmaxValues) || ismember(QnameLTA{q},QmaxValuesSpe)
                                maxValue = max(abs(value(:,:,a)));
                            else
                                maxValue = QmaxValuesMod{q};
                            end
                            value = value/maxValue;
                        end
                        
                        plot([startAPF:(endAPF-startAPF)/(Lt-1):endAPF], value(:,:,a), 'color', Qcolor, 'LineWidth', QLineWidth, 'LineStyle', QLineStyle);
                        
                        % Including x factor applied in the legend (3.4)
                        QfactorTag = '';
                        if strcmp(QplotRenorm,'renorm')% doesn't make sense to apply factor when plotRenormLTA = true (4.0)
                            QfactorTag = [' (' num2str(maxValue,2) ')']; % no need to put units in that case as they're all the same
                            ylim([-1,1]); % setting fixed limits
                        end
                        legendQname = [legendQname ; [legendAnimal{a} QfactorTag]]; %#ok<*AGROW> % 4.0
                        %--------------------------------------------------------------
                    end
                    
                    plotName = QuantityName;
                    tagProjectionFull = ['all.' QnameLTA{q} '_'  mapAnimal '.' uAnimal '_' tagProjectionTime  Qtags{t} '_' QplotType '_' QplotRenorm]; % using projectOR animal (4.0)

                else
                    %%% Plotting multiple values but only average values
                    %--------------------------------------------------------------
                    % plot, legend, color, units
                    % 4.0
                    if strcmp(QplotRenorm,'renorm')
                        if isempty(QmaxValues) || ismember(QnameLTA{q},QmaxValuesSpe)
                            maxValue = max(abs(value(:,:,1)));
                        else
                            maxValue = QmaxValuesMod{q};
                        end
                        value = value/maxValue;
                    end
                    
                    plot([startAPF:(endAPF-startAPF)/(Lt-1):endAPF], value(:,:,1), 'color', Qcolor, 'LineWidth', 2);
                    
                    % Including x factor applied in the legend (3.4)
                    QfactorTag = '';
                    if Qfactor(q) ~= 1 && strcmp(QplotRenorm,'raw') % doesn't make sense to apply factor when QplotRenorm = 'renorm' (4.0)
                        QfactorTag = [' x ' num2str(Qfactor(q)) ' (' Qunits ')'];
                    elseif strcmp(QplotRenorm,'renorm')
                        QfactorTag = [' (' num2str(maxValue,2) ' ' Qunits ')'];
                        ylim([-1,1]); % setting fixed limits
                    end
                    legendQname = [legendQname ; [plotName QfactorTag]]; % 3.5
                    %--------------------------------------------------------------
                    
                    plotName = mapAnimal;
                    tagProjectionFull = [mapAnimal '.' uAnimal '_' tagProjectionTime  Qtags{t} '_' QplotType '_' QplotRenorm]; % using projectOR animal (4.0)
                end
                
            else
                disp(['WARNING: quantity "' QuantityName '" was not found in LTA backup and was NOT plotted!'])
            end
        end
        
        legendQname = FormatLegend(legendQname); % replaces "dot" by ".u"
        
        %------------------------------------------------------------------
        % plot zeros for reference
        plot([startAPF:(endAPF-startAPF)/(Lt-1):endAPF],zeros(size([startAPF:(endAPF-startAPF)/(Lt-1):endAPF])),'-k','LineWidth',0.5); % mod 4.0
        hold off
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % title, legend, and other graphic plot tweak
        title([plotName ' ' Qtags{t}(2:end) ' ' num2str(uTimeWidth) 'h ' num2str(uTimeOverlap) ' olap | [' boxTag ']']);
        ylabel('Components');
        if length(QnameLTA) == 1 && strcmp(QplotRenorm,'raw')
           ylabel(['Components ( ' Qunits ' )']);
        end
        xlabel('Time ( hAPF )')
        legend(legendQname,'Location','eastoutside','FontSize',7)
        legend(legendQname,'Location','eastoutside')
        main_h = figure(1);
        set(main_h, 'color', 'white');
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % Print
        filename = [tagProjectionFull '_' boxTagFile ]; % 4.0
%         filename = ['PlotTime' Qtags{t} tagProjectionFull '_' boxTagFile ];
        if ~exist([OutputPathName filesep 'Plots_LTA'], 'dir')
            mkdir([OutputPathName filesep 'Plots_LTA']); 
        end
        if strcmp(imageExtension, '.svg')
            plot2svg([OutputPathName filesep 'Plots_LTA' filesep filename '.svg'],figure(1),'png')
        else
            print(printFormat, printResolution, [OutputPathName filesep 'Plots_LTA' filesep filename imageExtension]);
        end
        close
        %------------------------------------------------------------------
        fprintf(' Done\n')
    end
    
end % for each parts (dev, iso, ortho)
fprintf('\tDone\n')


%% History %%

% 13-25/06/2019: 4.0
% - many fixes and ergnonomy improvements (look for "4.0" tags)
% - fixed mistake in converting grid coordinates (grid center based) into
% grid indices.
