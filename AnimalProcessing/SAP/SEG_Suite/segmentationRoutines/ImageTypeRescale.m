function [newImage] = ImageTypeRescale(image,percentil,bitType,removeZeros)
%
% Function rescale the input image intensity distribution to a different
% a different bit depth and will saturate the extrema intensity value
% based on low and high percentil value
%
% Input
% rawImage    : the input image
% percentil   : [lowPer highPer] percentil value (default [0 100])
% bitDepth    : new image bit type to be cast in (default uint8)
% removeZeros : if true, will set zeros value to NaN (default false)
%
% example:
% newImage = BitDepthRescaler(rawImage,[lowPer highPer],bitDepth,zerosPix);
% newImage = BitDepthRescaler(rawImage,[1 99],'uint16',true);
%
% Stephane Rigaud
% 19/07/2016 - v1
%

% we make sure the input is good to be processed
doubleImage = double(image);
if removeZeros, doubleImage(doubleImage == 0) = NaN; end

% calculate distribution extrema to be removed
lowPerValue  = double(prctile(doubleImage(:), percentil(1)));
highPerValue = double(prctile(doubleImage(:), percentil(2)));

% cast image
doubleImage(doubleImage>highPerValue) = highPerValue;
doubleImage = doubleImage - lowPerValue;
newImage = doubleImage * (double(intmax(bitType))/highPerValue);
newImage = cast(newImage, bitType);

% if display
%     figure()
%     subplot(2,2,1)
%     imshow(uint16(image));
%     subplot(2,2,3)
%     hist(double(image));
%     subplot(2,2,2)
%     imshow(newImage);
%     subplot(2,2,4)
%     hist(double(newImage));
% end

end

%% Historic

% Todo: 
% - improve argument managment and flexibility
