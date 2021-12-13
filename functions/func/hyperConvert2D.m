function [Image2D] = hyperConvert2D(Image3D)
% Image3D is m*n*s
% Image2D is s*mn
if (ndims(Image3D) == 2)
    numBands = 1;
    [h, w] = size(Image3D);
else
    [h, w, numBands] = size(Image3D);
end
Image2D = reshape(Image3D, w*h, numBands).';
end

% aa(:,:,1)=[1,2;3,4;5,6];aa(:,:,2)=[1,2;3,4;5,6];bb=hyperConvert2D(aa)