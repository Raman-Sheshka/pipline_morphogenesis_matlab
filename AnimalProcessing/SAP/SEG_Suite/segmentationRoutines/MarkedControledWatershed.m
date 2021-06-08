function [W,bw,imgDist] = MarkedControledWatershed(I, S, sigma)

% I = imcomplement(I);
% bw = imbinarize(I,'global'); % read image


Ig1 = imgaussfilt(double(I),sigma);
Ig2 = imgaussfilt(double(I),sigma+1);

Ig = Ig1 - Ig2;
% bw = imbinarize(Ig,'global'); % read image
bw = imbinarize(Ig,'global');
% bw = imfill(bw,'holes');

bwSkel = ~bwmorph(~bw, 'Skel', Inf);
imgDist=-bwdist(~bwSkel);
imgDist=imimposemin(imgDist,S);

% imgDist=imimposemin(I,S);
W = ~watershed(imgDist);


end

