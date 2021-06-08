function I = Preprocessing(I)

background = imopen(I,strel('disk',15));
I = I - background;

% I = imadjust(I);
I = adapthisteq(I);

I = medfilt2(I);

I = imopen(I, strel('disk',1));

end

