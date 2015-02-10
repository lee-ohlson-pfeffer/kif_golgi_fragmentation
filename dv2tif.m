function [] = dv2tif(dir1, dirOut1, colorList1, wavelengthList1)

% dv2tif takes an input directory of projected dv images and creates
% seperate tif image files in the output directory based on the color and
% wavelengths given. These individual images can be used with CellProfiler
% for further processing.
%
% div2tif uses dv_data_extract_metadata which requires bio-formats
% (http://loci.wisc.edu/software/bio-formats)

%{
Here is an example of a colorList1, wavelengthList1 combination
colorList1 = {'GOL','RNA','KIF','NUC'};
wavelengthList1 = [676, 594, 523, 435];
%}

list1 = dir(dir1);
[list1Length1, ~] = size(list1);

mkdir(dirOut1);

parfor forA = 1:list1Length1
    
    name1 = list1(forA).name;
    
    if numel(name1)> 18 && isequal(name1(end-5:end-3),'PRJ')
        
        fileName1 = [dir1 list1(forA).name] 
        [imgIn1, metadata1] = dv_data_extract_metadata(fileName1);
        
        microscopeList1 = zeros(1,length(wavelengthList1));
        
        for forB = 1:length(microscopeList1)
            microscopeList1(forB) = metadata1.get(['Wavelength ' num2str(forB) ' (in nm)']);
        end
        
        colorList2 = cell(1,length(wavelengthList1));
        
        for forB = 1:length(colorList2)
            
            searchIndx1 = find(wavelengthList1==microscopeList1(forB),1);
            colorList2{forB} = colorList1{searchIndx1}; %#ok<PFBNS>
        end
   
        imgOut2 = uint16(imgIn1);
        
        for forB = 1:size(imgIn1,3)
            imwrite(imgOut2(:,:,forB), [dirOut1 list1(forA).name(1:end-14) colorList2{forB} '.tif']);
        end

    end
end