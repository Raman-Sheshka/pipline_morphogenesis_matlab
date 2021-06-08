function [ M ] = CreateRoiMask( I, W, R, T )

g = fspecial('gaussian',[R R],round(R/3));
I = im2double(I);
I = imadjust(I,[nanmin(I(:)) nanmax(I(:))], []);

% variance filter and threshold
I = stdfilt(I, ones(W));
I = imfilter(I, g, 'replicate');
level = graythresh(I);
M = im2bw(I,level/4);

% morphomath cleaning
M = imdilate(M, strel('disk',10));
M = imerode(M, strel('disk',30));
% M = imclose(M, strel('disk',10));
% M = imopen(M, strel('disk',20));
M = imfill(M, 'holes');

% small region cleaning
CC = bwconncomp(M);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,~] = max(numPixels);
M = bwareaopen(M, biggest-1);
CC = bwconncomp(~M);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,~] = max(numPixels);
M = ~bwareaopen(~M, round(biggest/100));

% M = imerode(M, strel('disk',20));
% g = fspecial('gaussian',75,round(75/3));
% M = imfilter(double(M), g, 'replicate');
% level = graythresh(M);
% M = im2bw(M,level);

end

