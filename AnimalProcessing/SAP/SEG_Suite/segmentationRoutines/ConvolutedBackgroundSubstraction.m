function [output]=ConvolutedBackgroundSubstraction(Img,radius)

if isa(Img,'uint8')
    type = 'uint8';
elseif isa(Img,'uint16')
    type = 'uint16';
else
    type = 'uint16';
end

diameter = round(radius.*2);
fun = @(x) MedianSubstraction(x(:));
output = nlfilter(Img,[diameter diameter],fun);

output = cast(output,type);

end

function [output]=MedianSubstraction(X)
index = round(length(X)./2);
output = single(X(index)) - median(single(X));
end