% Local Time Analysis
%
% version 3
% Stéphane Rigaud

% Plot the evolution of quantities (tensor, vector, scalar) over time in a
% specific area of the animal.
%
% if the GUI boolean is true, display the animal grid, with landmark
% otherwise, use the provided coordinates. If no coordinates provided, we
% plot the largest possible area.
%
% todo: futur update, add the std values in the plot when relevant.


% If GUI, use graphic interface to define zone.
% Only one zone will be treated this way, for multiple zone, use Qareas
% if GUI
%     Qareas = cell(1);
%     %%% BG image
%
%     BACKUP.overlap = 0;
%     size_image_x = BACKUP.xywh(3)*(BACKUP.size(2)*BACKUP.overlap+1);
%     size_image_y = BACKUP.xywh(4)*(BACKUP.size(1)*BACKUP.overlap+1);
%     if BACKUP.overlap == 0
%         size_image_x = BACKUP.xywh(3)*(BACKUP.size(2)+1);
%         size_image_y = BACKUP.xywh(4)*(BACKUP.size(1)+1);
%         BACKUP.overlap=1;
%     end
%     Qimage = imresize(Qimage,'Method','bicubic','OutputSize',[size_image_y size_image_x]);
%     %%% zone selection
%     waitSelectionLoop = true;
%     while(waitSelectionLoop)
%         imshow(Qimage, 'Border','tight');
%         % Get XY coordinate of the zone
%         hold on
%         [xi,yi,buttemp]= ginput(2);
%         xi=[xi(1) xi(2) xi(2) xi(1)];
%         yi=[yi(1) yi(1) yi(2) yi(2)];
%         X = round(xi./ (BACKUP.xywh(3)*BACKUP.overlap));
%         X(X==0) = 1;
%         X(X>BACKUP.size(2)) = BACKUP.size(2);
%         Y = round(yi./ (BACKUP.xywh(4)*BACKUP.overlap));
%         Y(Y==0) = 1;
%         Y(Y>BACKUP.size(1)) = BACKUP.size(1);
%         nX = [X X(1)] .* (BACKUP.xywh(3)*BACKUP.overlap);
%         nY = [Y Y(1)] .* (BACKUP.xywh(3)*BACKUP.overlap);
%         plot(nX,nY,'-y')
%         hold off
%         % Check with user
%         button = questdlg('Is this the region you want to process?','Selection Validation','Yes','No','Yes');
%         if strcmp(button,'Yes')
%             waitSelectionLoop = false;
%             CropRegionGrid = [X(1),X(2);Y(2),Y(3)];
%             if saveCropRegion
%                 filename = ['Contributions_Crop_' Animal '_(' num2str(CropRegionGrid(1,1)) ',' num2str(CropRegionGrid(2,1)) ')(' num2str(CropRegionGrid(1,2)) ',' num2str(CropRegionGrid(2,2)) ')'];
%                 if ~exist([OutputPathName filesep 'Plots'],'dir'); mkdir([OutputPathName filesep 'Plots']); end;
%                 print('-dpng', '-r100', [OutputPathName filesep 'Plots' filesep filename '.png']);
%             end
%             close;
%         else
%             close;
%         end
%     end
%     % Save zone
%     Qareas{1} = CropRegionGrid;
% end

[Ly, Lx, Lz, Lt] = size(eval(['BACKUP.' Qname{1} tagTensor{1}]));
if isempty(Qareas)
    disp(['WARNING: no region specified, we take the largest region possible -> ' num2str([1 Lx 1 Ly])]);
    Qareas = {[1 Lx ; 1 Ly]};
end

%------------------------------------------------------------------
% Absisse Initialisation
startFrame = (BACKUP.FrameArray(1,1) + BACKUP.FrameArray(1,2))   ./ 2;
endFrame   = (BACKUP.FrameArray(Lt,1) + BACKUP.FrameArray(Lt,2)) ./ 2;
startAPF   = round( frame2time(startFrame, BACKUP.TimeArray{1,1}, BACKUP.FrameArray(1,1), delta_t, 'dec') );
endAPF     = round( frame2time(endFrame, BACKUP.TimeArray{1,1}, BACKUP.FrameArray(1,1), delta_t, 'dec') );
%------------------------------------------------------------------

for t = 1:length(tagTensor)
    for z = 1:length(Qareas)
        
        [Ly, Lx, Lz, Lt] = size(eval(['BACKUP.' Qname{1} tagTensor{1}]));
        Total = zeros(1,Lt);
        Qmean = cell(1,size(Qname,1));
        legendQname = cell(size(Qname));
        CropRegionGrid = Qareas{z};
        
        %% Quantity extraction --------------------------------------------
        for q = 1:length(Qname)  % for each quantities
                           
            %--------------------------------------------------------------
            % Get information the data
            QuantityName = [ Qname{q} tagTensor{t}];
            [Q1,Q2] = QuantityNameSpliter(QuantityName);
            Pname = GetPname(Q1);
            uPname = '';
            if ~isempty(Q2)
                uPname = GetPname(Q2);
            end
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % weight calculation
            Weight  = eval(['BACKUP.AreaRatios_' Pname '(CropRegionGrid(2,1):CropRegionGrid(2,2),CropRegionGrid(1,1):CropRegionGrid(1,2),:,:);']);
            if ~isempty(Q2)
                uAR = eval(['BACKUP.AreaRatios_' uPname '(CropRegionGrid(2,1):CropRegionGrid(2,2),CropRegionGrid(1,1):CropRegionGrid(1,2),:,:);']);
                Weight = min(Weight, uAR);
            end
            if strcmp(Pname, 'TA') || strcmp(uPname, 'TA')
                RConds = BACKUP.RConds(CropRegionGrid(2,1):CropRegionGrid(2,2),CropRegionGrid(1,1):CropRegionGrid(1,2),:,:);
                Weight = Weight .* RConds;
                Weight(isnan(Weight)) = 0;
            end
            Weight = Weight .^ 2;
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % Calculate qunatity to plot
            Quantity = eval(['BACKUP.' QuantityName ';']);
            if ndims(Quantity) == 3
                [qy, qx, ~, qt] = size(Quantity);
                Quantity = reshape(Quantity,[Ly, Lx, 1, Lt]);
            end
            QuantityCrop = Quantity(CropRegionGrid(2,1):CropRegionGrid(2,2),CropRegionGrid(1,1):CropRegionGrid(1,2),:,:);
            if size(QuantityCrop,3) ~= 1 % not a scalar
               QuantityCrop = TensorNormMap(QuantityCrop);
            end
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % Calculate pondarate mean of the crop zone
            regionY = length(CropRegionGrid(2,1):CropRegionGrid(2,2));
            regionX = length(CropRegionGrid(1,1):CropRegionGrid(1,2));
            resWeight = reshape(Weight,[regionY regionX 1 Lt]);
            weightedQuantityCrop = QuantityCrop .* resWeight;
            weightedMeanCrop = nansum(nansum(weightedQuantityCrop,1),2) ./ sum(sum(resWeight,1),2);
            weightedMeanCrop = reshape(weightedMeanCrop,[1 Lt]);
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % Normalisation, <||Q||>
            if strcmp(Normalisation,'animal') % over the animal
                QuantityNorm = Quantity;
                if size(QuantityNorm,3) ~= 1 % not a scalar
                    QuantityNorm = TensorNormMap(QuantityNorm);
                end
                Weight = eval(['BACKUP.AreaRatios_' Pname ';']);
                if ~isempty(Q2)
                    uAR = eval(['BACKUP.AreaRatios_' uPname ';']);
                    Weight = min(Weight, uAR);
                end
                if strcmp(Pname, 'TA') || strcmp(uPname, 'TA')
                    Weight = Weight .* BACKUP.RConds;
                    Weight(isnan(Weight)) = 0;
                end
                Weight = Weight .^ 2;
                weightedQuantityNorm = QuantityNorm .* reshape(Weight,[Ly Lx 1 Lt]);
                weightedMeanNorm = nansum(weightedQuantityNorm(:)) ./ sum(Weight(:));
                weightedMeanCrop = weightedMeanCrop ./ weightedMeanNorm;
            elseif strcmp(Normalisation,'region') % over the region
                weightedMeanNorm = mean(weightedMeanCrop(:));
                weightedMeanCrop = weightedMeanCrop ./ weightedMeanNorm;
            end
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % save
            Qmean{q} = weightedMeanCrop;
            Total = Total + weightedMeanCrop;
            %--------------------------------------------------------------
            
        end  % end for each quantities
        %------------------------------------------------------------------
        
        %% Plot -----------------------------------------------------------
        figure();
        if ~isempty(Qrange)
            axis([ startAPF-1 endAPF+1 Qrange{t}]);
        end
        hold on
        for q = 1:length(Qname)
            QuantityName = [ Qname{q} tagTensor{t}];
            
            %--------------------------------------------------------------
            % manage complex quantity name
            [Q1,Q2] = QuantityNameSpliter(QuantityName);
            plotName = QuantityName;
            [Pname,idx] = GetPname(Q1);
            uPname = '';
            if ~isempty(Q2)
                uPname = GetPname(Q2);
            end
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % get value and mean value
            value = Qmean{q};
            if length(value) == 1
                startAPF   = round( Time_str2dec(TIME.multiTimeStart) );
                endAPF     = round( Time_str2dec(TIME.multiTimeStop) );
                Lt = 2;
                value = Qmean{q} .* ones( size( [startAPF:(endAPF-startAPF)/(Lt-1):endAPF] ) );
            end
            %--------------------------------------------------------------
            
            %--------------------------------------------------------------
            % plot, legend, color, units
            Qcolor = eval(['allColors_' Pname '{idx}']);
            Qunits = eval(['allUnits_'  Pname '{idx}']);
            plot( [startAPF:(endAPF-startAPF)/(Lt-1):endAPF], value, 'color', Qcolor, 'LineWidth', 2);
            legendQname{q} = plotName;
            %--------------------------------------------------------------
        end
        
        %------------------------------------------------------------------
        % plot zeros for reference
        plot([startAPF:(endAPF-startAPF)/(Lt-1):endAPF],zeros(size([startAPF:(endAPF-startAPF)/(Lt-1):endAPF])),'-.k','LineWidth',2);
        hold off
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % title, legend, and other graphic plot tweak
        title(['Plot ' Animal ' ' tagTensor{t}(2:end) ' ' num2str(uTimeWidth) 'h ' num2str(uTimeOverlap) ' olap | [(' num2str(CropRegionGrid(1,1)) ',' num2str(CropRegionGrid(2,1)) ');(' num2str(CropRegionGrid(1,2)) ',' num2str(CropRegionGrid(2,2)) ')]']);
        ylabel(['Processes ( ' Qunits ' )'])
        xlabel('Time ( hAPF )')
        legend(legendQname,'Location','eastoutside')
        main_h = figure(1);
        set(main_h, 'color', 'white');
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % Print 
        filename = ['PlotTime' tagTensor{t} tagProjection '_' Animal '_(' num2str(CropRegionGrid(1,1)) ',' num2str(CropRegionGrid(2,1)) ')(' num2str(CropRegionGrid(1,2)) ',' num2str(CropRegionGrid(2,2)) ')'];
        if ~exist([OutputPath filesep 'Plots_LTA'], 'dir'); mkdir([OutputPath filesep 'Plots_LTA']); end;
        if strcmp(imageExtension, '.svg')
            plot2svg([OutputPath filesep 'Plots_LTA' filesep filename '.svg'],figure(1),'png')
        else
            print(printFormat, printResolution, [OutputPath filesep 'Plots_LTA' filesep filename imageExtension]);
        end
        close
        %------------------------------------------------------------------
        
    end % for each zones
end % for each parts (dev, iso, ortho)

% 26/05/206: 3.1
% - weight verification to remove all NaN coming from RConds

% 13/03/2016: 3.0