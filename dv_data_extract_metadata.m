function [imgOut1, metadataOut1] = dv_data_extract_metadata(fileName1)

% Helper function for opening dv files and extracting metadata. Uses
% bio-formats (http://loci.wisc.edu/software/bio-formats)

imgIn1 = bfopen(fileName1);
imgIn2 = imgIn1{1,1};
lengthImgIn2 = size(imgIn2,1);
sizeImgIn2Slice1 = size(imgIn2{1,1});
imgOut1 = zeros(sizeImgIn2Slice1(1), sizeImgIn2Slice1(2), lengthImgIn2);
for forA = 1:lengthImgIn2
    imgOut1(:,:,forA) = imgIn2{forA,1};
end

metadataOut1 = imgIn1{2};