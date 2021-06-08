function DBAbackup = MakeDBAbackup(LoadedBackupsIn, DISPLAY, gridOverlap)
%
% DBAbackup = MakeDBAbackup(LoadedBackupsIn, DISPLAY, gridOverlap)
%
% Version 1.0
% Boris Guirao


%% Code %%

% Get the grid size
s1 = LoadedBackupsIn{1}.Size; % (H1,W1)
s2 = LoadedBackupsIn{2}.Size; % (H2,W2)

oIJ1 = LoadedBackupsIn{1}.originBoxIJ;
oIJ2 = LoadedBackupsIn{2}.originBoxIJ;

oIJc = max([oIJ1;oIJ2],[],1); %  originBoxIJ in common grid
sc = oIJc + max([s1-oIJ1; s2-oIJ2],[],1); % common grid size (Hc,Wc)

%% Determining startLocIJ, endLocIJ %%

startLocIJ = 1 + oIJc - [oIJ1 ; oIJ2];
endLocIJ = startLocIJ + [s1 ; s2] - 1;

%% Determining image size

% Minimal image size
sizeImageXmin = ceil(DISPLAY.boxSize(1) * (1 + (sc(2)-1) * (1 - gridOverlap)));
sizeImageYmin = ceil(DISPLAY.boxSize(2) * (1 + (sc(1)-1) * (1 - gridOverlap)));

% Adding at least one full box size along X and Y on BOTH sides of images (2.7)
leeway = 0.05;
sizeImageXext = ceil(min(sizeImageXmin * (1+leeway), sizeImageXmin + 2 * DISPLAY.boxSize(1)));
sizeImageYext = ceil(min(sizeImageYmin * (1+leeway), sizeImageYmin + 2 * DISPLAY.boxSize(2)));
% NB: Adding leeway so as NOT to start grid at [0,0] but further from image border and still fit it in image

% setting final leeway:
leewayPixX = sizeImageXext - sizeImageXmin; % image extra pixels along X
leewayPixY = sizeImageYext - sizeImageYmin; % image extra pixels along Y

sizeImageXreg = sizeImageXext;
sizeImageYreg = sizeImageYext;

% Setting new origin
% % ACTUAL new origin in image:
newOxReg = DISPLAY.boxSize(1) * (oIJc(2)-1) * (1 - gridOverlap) + leewayPixX/2;
newOyReg = DISPLAY.boxSize(2) * (oIJc(1)-1) * (1 - gridOverlap) + leewayPixY/2;
% meanYmidReg = newOyReg - meanDeltaYmid; % 2.5
% NB: subtracting "PLOT.boxSize" because box origin is normally at ULC of box => subtraction of one full box width/height

% Build global common REGULAR grid
% since using gridSize to keep previous size, need to determine grid new ULC in extended image (2.7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRIDreg = MakeGrid([sizeImageYreg sizeImageXreg], DISPLAY.boxSize, [leewayPixX/2 leewayPixY/2], sc , [], [], gridOverlap);
% NB: at this stage, array "Coordinates" is NOT up to date as temporary
% grid origin was taken at top left corner of grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updating grid xywh by actual origin [newOxReg newOyReg]:
GRIDreg.xywh = [newOxReg newOyReg GRIDreg.xywh(3) GRIDreg.xywh(4)];
% Updating array "Coordinates" to set compartment containing orgin to [0,0] (3.4):
[gridCoordinates, originBoxIJcheck] = gridULCs2gridCoordinates(GRIDreg.ULCs, round([newOxReg newOyReg]));

if norm(originBoxIJcheck - oIJc) > 0
    warndlg('Difference of "originBoxIJ"!')
end

% Filling REG structure
REG.xywh = GRIDreg.xywh;
REG.Centroids = GRIDreg.Centroids;
REG.SizeImageX = sizeImageXreg;
REG.SizeImageY = sizeImageYreg;


%% Filling "DBAbackup" %%

fields2skip = {'Coordinates', 'ULCs', 'Color','FrameArray','TimeArray','LineWidth','Size','fullImage','originBoxIJ','Overlap'};

fieldList = fieldnames(LoadedBackupsIn{1}); % Get the list of field to process from one of the backup
for f = 1:length(fieldList)
    
    fthFieldIn1 = LoadedBackupsIn{1}.(fieldList{f});
    fthFieldIn2 = LoadedBackupsIn{2}.(fieldList{f});
    
    if ~ismember(fieldList{f},fields2skip) && ~isstruct(fthFieldIn1)
        
        fthFieldOut1 = zeros(sc(1),sc(2),size(fthFieldIn1,3),size(fthFieldIn1,4),size(fthFieldIn1,5));
        fthFieldOut1(startLocIJ(1,1):endLocIJ(1,1), startLocIJ(1,2):endLocIJ(1,2),:,:,:) = fthFieldIn1;
        fthFieldOut1 = fthFieldOut1(:,:,:,:,1); % only keeps average value (not each animal value)
        
        fthFieldOut2 = zeros(sc(1),sc(2),size(fthFieldIn2,3),size(fthFieldIn2,4),size(fthFieldIn2,5));
        fthFieldOut2(startLocIJ(2,1):endLocIJ(2,1), startLocIJ(2,2):endLocIJ(2,2),:,:,:) = fthFieldIn2;
        fthFieldOut2 = fthFieldOut2(:,:,:,:,1); % only keeps average value (not each animal value)
        
        if strcmp(fieldList{f},'AreaRatios') % min(A.ar, B.ar)
            disp(['min : ' fieldList{f}]);
            ARfield = min(fthFieldOut1, fthFieldOut2);
            DBAbackup.(fieldList{f}) = ARfield;
            
        elseif  ~isempty(strfind(fieldList{f},'_Smap')) % if Significance map, logical AND(A.map, B.map)
            disp(['logical OR : ' fieldList{f}]);
            fthFieldSmap = or(fthFieldOut1, fthFieldOut2);
            DBAbackup.(fieldList{f}) = fthFieldSmap;
        else
            fthDeltaFields = fthFieldOut1 - fthFieldOut2;
            DBAbackup.(fieldList{f}) = fthDeltaFields;
        end
        
    elseif ~ismember(fieldList{f},fields2skip) && isstruct(fthFieldIn1) % REG & FULL structures
        
        if strcmp(fieldList{f},'REG')
            
            % updating macro coordinates in new common image
            macro1c = fthFieldIn1.Macrocaetes - repmat(oIJ1',1,8) + repmat(oIJc',1,8);
            macro2c = fthFieldIn2.Macrocaetes - repmat(oIJ2',1,8) + repmat(oIJc',1,8);
            % taking average
            REG.Macrocaetes = 0.5*(macro1c + macro2c);
            
            % updating midline coordinates in new common image
            yMid1c = fthFieldIn1.yMid - oIJ1(1) + oIJc(1);
            yMid2c = fthFieldIn2.yMid - oIJ2(1) + oIJc(1);
            % taking average
            REG.yMid = 0.5*(yMid1c + yMid2c);
            
            DBAbackup.REG = REG;
            
        elseif strcmp(fieldList{f},'FULL')
            
            FULL = struct();
            DBAbackup.FULL = FULL;
        end
        
    else
        DBAbackup.(fieldList{f}) = fthFieldIn1; % copy FIRST value of fields2skip as it is
    end
end

%% Overwritting %%

%common size:
DBAbackup.Size = sc;
% originBoxIJ
DBAbackup.originBoxIJ = oIJc;
% ULCs, Coordinates
DBAbackup.Coordinates = gridCoordinates;
DBAbackup.ULCs = GRIDreg.ULCs;



%% History %%

% 03/06/20: creation

