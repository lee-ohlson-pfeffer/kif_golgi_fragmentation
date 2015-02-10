function [cel_stat,cel_legn]= postcp_sirna(fileoutmat)

% postcp_sirna takes the mat output from cellprofiler
% (http://www.cellprofiler.org/) and extracts relevant feature properties.
% This function works for siRNA analysis. The nocodoazole analysis uses the
% same function lacking the rna intensity metrics.

load(fileoutmat);

%{
The following parameters were used for calculating area thresholds
%}

small_area_thresh = 15;
large_distance_thresh = 2.5;
object_thresh_percent = .1;
object_thresh_thresh = 100;

DirOut1 = [handles.Current.DefaultOutputDirectory '\'];

img_names = handles.Measurements.Image.FileName_NUC';
img_boxs = str2double(cellfun((@(x) x(1)),img_names,'UniformOutput',false));
img_slid = str2double(cellfun((@(x) x(3:4)),img_names,'UniformOutput',false));
img_slip = str2double(cellfun((@(x) x(6)),img_names,'UniformOutput',false));
img_site = str2double(cellfun((@(x) x(8:10)),img_names,'UniformOutput',false));
img_prop = [img_boxs img_slid img_slip img_site];

nuc_numb = cell2mat(handles.Measurements.Image.Count_Nuclei)';
tot_num = sum(nuc_numb);
cel_stat = double(zeros(tot_num,4));
count1 = 1;
for for1 = 1:numel(nuc_numb)
    cel_stat(count1:count1+nuc_numb(for1)-1,:) = repmat(img_prop(for1,:),nuc_numb(for1),1);
    count1 = count1+nuc_numb(for1);
end

file_name = handles.Measurements.Image.FileName_GOL;

cel_numb = [nuc_numb ones(size(nuc_numb))];
cel_numb = mat2cell(cel_numb,ones(size(nuc_numb)))';

cgl_area = cell(numel(file_name),1);
cgl_conv = cell(numel(file_name),1);
cgl_avgI = cell(numel(file_name),1);
cgl_intI = cell(numel(file_name),1);
cgl_eqdi = cell(numel(file_name),1);
cgl_soli = cell(numel(file_name),1);
cgl_mrad = cell(numel(file_name),1);
cgl_vrad = cell(numel(file_name),1);
cgl_miax = cell(numel(file_name),1);
cgl_maax = cell(numel(file_name),1);
cgl_orie = cell(numel(file_name),1);
cgl_mdis = cell(numel(file_name),1);
cgl_vdis = cell(numel(file_name),1);

bad_mgls = cell(numel(file_name),1);

for for1 = 1:numel(file_name)
    load([DirOut1 file_name{for1}(1:end-7) 'GOL.mat'])
    gol_obj = Image;
    load([DirOut1 file_name{for1}(1:end-7) 'NUC.mat'])
    gol_img = Image;
    load([DirOut1 file_name{for1}(1:end-7) 'RNA.mat'])
    gol_msk = Image;
    
    temp_nuc_num = nuc_numb(for1);
        
    gol_obj_mod = zeros(size(gol_obj));
    
    for for2 = 1:temp_nuc_num
        temp_gol_obj = (gol_obj==for2);
        if sum(temp_gol_obj(:))>0
            toggle1 = 1;
            while toggle1 == 1
                temp_stat = regionprops(double(temp_gol_obj), gol_img, 'WeightedCentroid');
                temp_stat_part = regionprops(temp_gol_obj, gol_img, 'PixelIdxList', 'Area', 'WeightedCentroid');
                temp_cent = temp_stat.WeightedCentroid;
                temp_list = cat(1,temp_stat_part.WeightedCentroid);
                temp_area = cat(1,temp_stat_part.Area);
                
                temp_cent = repmat(temp_cent,size(temp_list,1),1);
                temp_diff = temp_list - temp_cent;
                temp_diff = temp_diff.^2;
                temp_diff = sum(temp_diff,2);
                temp_diff = (temp_diff.^0.5);
                
                temp_vari = var(temp_diff,temp_area);
                temp_stde = temp_vari.^.5;
                temp_mean = mean(temp_diff);
                
                temp_mdif = (temp_diff - temp_mean)./temp_stde;
                [temp_maxx,temp_indx] = max(temp_mdif);
                
                if temp_maxx>large_distance_thresh && temp_area(temp_indx)<small_area_thresh
                    temp_gol_obj(temp_stat_part(temp_indx).PixelIdxList) = 0;
                else
                    toggle1 = 0;
                end
            end
            gol_obj_mod(temp_gol_obj) = for2;
        end
    end
    
    temp_min_obj = gol_msk((gol_msk>0)&(gol_obj_mod>0));
    temp_min_obj = unique(temp_min_obj);
    
    temp_stat = regionprops(gol_obj_mod,gol_img, 'all');

    im_gol_mrad = NaN(temp_nuc_num,1);
    im_gol_vrad = NaN(temp_nuc_num,1);
    im_gol_mdis = NaN(temp_nuc_num,1);
    im_gol_vdis = NaN(temp_nuc_num,1);
    im_gol_area = NaN(temp_nuc_num,1);
    im_gol_conv = NaN(temp_nuc_num,1);
    im_gol_avgI = NaN(temp_nuc_num,1);
    im_gol_intI = NaN(temp_nuc_num,1);
    im_gol_eqdi = NaN(temp_nuc_num,1);
    im_gol_soli = NaN(temp_nuc_num,1);
    im_gol_miax = NaN(temp_nuc_num,1);
    im_gol_maax = NaN(temp_nuc_num,1);
    im_gol_orie = NaN(temp_nuc_num,1);
    
    im_gol_area(1:length(temp_stat)) = cat(1,temp_stat.Area);
    im_gol_conv(1:length(temp_stat)) = cat(1,temp_stat.ConvexArea);
    im_gol_avgI(1:length(temp_stat)) = cat(1,temp_stat.MeanIntensity);
    im_gol_intI(1:length(temp_stat)) = cat(1,temp_stat.MeanIntensity).*cat(1,temp_stat.Area);
    im_gol_eqdi(1:length(temp_stat)) = cat(1,temp_stat.EquivDiameter);
    im_gol_soli(1:length(temp_stat)) = cat(1,temp_stat.Solidity);
    im_gol_miax(1:length(temp_stat)) = cat(1,temp_stat.MinorAxisLength);
    im_gol_maax(1:length(temp_stat)) = cat(1,temp_stat.MajorAxisLength);
    im_gol_orie(1:length(temp_stat)) = cat(1,temp_stat.Orientation);
    
    parfor for2 = 1:numel(temp_stat)
        temp_list = temp_stat(for2).PixelList;
        temp_valu = temp_stat(for2).PixelValues;
        temp_cent = temp_stat(for2).WeightedCentroid;
        temp_intI = sum(temp_valu);
        
        if numel(temp_list)>2
            pcaCoeff = princomp(temp_list);
        else
            pcaCoeff = eye(2);
        end
        
        temp_rlis = temp_list;
        temp_rlis(:,1) = (temp_list(:,1)-temp_cent(1,1)).*pcaCoeff(1)+(temp_list(:,2)-temp_cent(1,2)).*pcaCoeff(2);
        temp_rlis(:,2) = (temp_list(:,1)-temp_cent(1,1)).*pcaCoeff(3)+(temp_list(:,2)-temp_cent(1,2)).*pcaCoeff(4);
        temp_mdis = sum(abs(temp_rlis(:,2)).*temp_valu)./temp_intI;
        temp_vdis = var(temp_rlis(:,2),temp_valu);
        
        im_gol_mdis(for2) = temp_mdis;
        im_gol_vdis(for2) = temp_vdis;
        
        temp_cent = repmat(temp_cent,numel(temp_valu),1);
        temp_diff = temp_list - temp_cent;
        temp_diff = temp_diff.^2;
        temp_diff = sum(temp_diff,2);
        temp_diff = (temp_diff.^0.5);
        
        temp_mrad = sum(temp_diff.*temp_valu)./temp_intI;
        temp_vrad = sum((temp_diff-temp_mrad).^2.*temp_valu)./temp_intI; % = var(temp_diff,temp_valu);
        
        im_gol_mrad(for2) = temp_mrad;
        im_gol_vrad(for2) = temp_vrad;
        
    end
    
    cgl_area{for1} = im_gol_area;
    cgl_conv{for1} = im_gol_conv;
    cgl_avgI{for1} = im_gol_avgI;
    cgl_intI{for1} = im_gol_intI;
    cgl_eqdi{for1} = im_gol_eqdi;
    cgl_soli{for1} = im_gol_soli;
    cgl_mrad{for1} = im_gol_mrad;
    cgl_vrad{for1} = im_gol_vrad;
    cgl_miax{for1} = im_gol_miax;
    cgl_maax{for1} = im_gol_maax;
    cgl_orie{for1} = im_gol_orie;
    cgl_mdis{for1} = im_gol_mdis;
    cgl_vdis{for1} = im_gol_vdis;
    
    bad_mgls{for1} = temp_min_obj;
end

nuc_labl = double(cell2mat(handles.Measurements.Nuclei.Number_Object_Number'));
nuc_area = double(cell2mat(handles.Measurements.Nuclei.AreaShape_Area'));
nuc_comp = cell2mat(handles.Measurements.Nuclei.AreaShape_Compactness');
nuc_soli = cell2mat(handles.Measurements.Nuclei.AreaShape_Solidity');
nuc_form = cell2mat(handles.Measurements.Nuclei.AreaShape_FormFactor');

nuc_perc = cell2mat(handles.Measurements.Nuclei.Neighbors_PercentTouching_Adjacent');
nuc_locy = cell2mat(handles.Measurements.Nuclei.Location_Center_Y');
nuc_locx = cell2mat(handles.Measurements.Nuclei.Location_Center_X');
nuc_miax = cell2mat(handles.Measurements.Nuclei.AreaShape_MinorAxisLength');
nuc_maax = cell2mat(handles.Measurements.Nuclei.AreaShape_MajorAxisLength');
nuc_mifd = cell2mat(handles.Measurements.Nuclei.AreaShape_MinFeretDiameter');
nuc_mafd = cell2mat(handles.Measurements.Nuclei.AreaShape_MaxFeretDiameter');

cel_area = double(cell2mat(handles.Measurements.Cells.AreaShape_Area'));
cel_kif_intI = cell2mat(handles.Measurements.Cells.Intensity_IntegratedIntensity_CropKIF');
cel_kif_avgI = cell2mat(handles.Measurements.Cells.Intensity_MeanIntensity_CropKIF');
cel_kif_medI = cell2mat(handles.Measurements.Cells.Intensity_MedianIntensity_CropKIF');

cel_kif_maxI = cell2mat(handles.Measurements.Cells.Intensity_MaxIntensity_CropKIF');
cel_kif_minI = cell2mat(handles.Measurements.Cells.Intensity_MinIntensity_CropKIF');
cel_kif_loQI = cell2mat(handles.Measurements.Cells.Intensity_LowerQuartileIntensity_CropKIF');
cel_kif_upQI = cell2mat(handles.Measurements.Cells.Intensity_UpperQuartileIntensity_CropKIF');
cel_mafd = cell2mat(handles.Measurements.Cells.AreaShape_MaxFeretDiameter');
cel_mifd = cell2mat(handles.Measurements.Cells.AreaShape_MinFeretDiameter');
cel_miax = cell2mat(handles.Measurements.Cells.AreaShape_MinorAxisLength');
cel_maax = cell2mat(handles.Measurements.Cells.AreaShape_MajorAxisLength');
cel_orie = cell2mat(handles.Measurements.Cells.AreaShape_Orientation');
cel_rna_loQI = cell2mat(handles.Measurements.Cells.Intensity_LowerQuartileIntensity_CropRNA');
cel_rna_minI = cell2mat(handles.Measurements.Cells.Intensity_MinIntensity_CropRNA');

bad_mgls = bad_mgls';
mgl_cel_pare = handles.Measurements.MaskedGolgi.Parent_Cells;
mgl_area = handles.Measurements.MaskedGolgi.AreaShape_Area;

mgl_cel_pare = cellfun(@(x,y) x(y), mgl_cel_pare, bad_mgls, 'UniformOutput',false);
mgl_area = cellfun(@(x,y) x(y), mgl_area, bad_mgls, 'UniformOutput',false);

cel_gol_numb = cellfun(@(x,z) accumarray(x,1,z), mgl_cel_pare, cel_numb, 'UniformOutput',false);
cel_gol_numb = cell2mat(cel_gol_numb');

cel_mgl_avgA = cellfun(@(x,y,z) accumarray(x,y,z,@mean), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);
cel_mgl_avgA = cell2mat(cel_mgl_avgA');
cel_mgl_avgA(cel_mgl_avgA==0) = NaN;

cel_mgl_maxA = cellfun(@(x,y,z) accumarray(x,y,z,@max), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);
cel_mgl_maxA = double(cell2mat(cel_mgl_maxA'));
cel_mgl_maxA(isnan(cel_mgl_avgA)) = NaN;

hist_bins = 10.^[0:0.4:4.4]; %#ok<NBRAK>
mgl_dis_area = cellfun(@(x,y,z) accumarray(x,y,z, @(w) {(histc(double(w),hist_bins,1)')./sum(histc(double(w),hist_bins,1))},{NaN(size(hist_bins))}), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);
mgl_dis_area = cat(1,mgl_dis_area{:});
mgl_dis_area = cat(1,mgl_dis_area{:});

mgl_100_area = cellfun(@(x,y,z) accumarray(x,y,z, @(w) sum(w(w>100))./sum(w),NaN), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);
mgl_150_area = cellfun(@(x,y,z) accumarray(x,y,z, @(w) sum(w(w>150))./sum(w),NaN), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);
mgl_225_area = cellfun(@(x,y,z) accumarray(x,y,z, @(w) sum(w(w>225))./sum(w),NaN), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);
mgl_338_area = cellfun(@(x,y,z) accumarray(x,y,z, @(w) sum(w(w>337.5))./sum(w),NaN), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);

mgl_per_area = cellfun(@(x,y,z) accumarray(x,y,z, @(w) sum(w(w>max(sum(w)*object_thresh_percent,object_thresh_thresh)))./sum(w),NaN), mgl_cel_pare, mgl_area, cel_numb, 'UniformOutput',false);

mgl_100_area = cell2mat(mgl_100_area');
mgl_150_area = cell2mat(mgl_150_area');
mgl_225_area = cell2mat(mgl_225_area');
mgl_338_area = cell2mat(mgl_338_area');
mgl_per_area = cell2mat(mgl_per_area');

mgl_intI = cell2mat(handles.Measurements.Cells.Mean_MaskedGolgi_Intensity_IntegratedIntensity_CropGOL');
mgl_maxI = cell2mat(handles.Measurements.Cells.Mean_MaskedGolgi_Intensity_MaxIntensity_CropGOL');
mgl_avgI = cell2mat(handles.Measurements.Cells.Mean_MaskedGolgi_Intensity_MeanIntensity_CropGOL');

cel_size = nuc_numb;
cel_size = mat2cell(cel_size,ones(size(nuc_numb)))';

cgl_comp = (handles.Measurements.CellGolgi.AreaShape_Compactness);
cgl_comp = cellfun(@(x,y) padarray(x,y-length(x),NaN,'post'), cgl_comp, cel_size, 'UniformOutput',false);
cgl_comp = double(cell2mat(cgl_comp'));

cgl_form = (handles.Measurements.CellGolgi.AreaShape_FormFactor);
cgl_form = cellfun(@(x,y) padarray(x,y-length(x),NaN,'post'), cgl_form, cel_size, 'UniformOutput',false);
cgl_form = double(cell2mat(cgl_form'));

cgl_mifd = handles.Measurements.CellGolgi.AreaShape_MinFeretDiameter;
cgl_mifd = cellfun(@(x,y) padarray(x,y-length(x),NaN,'post'), cgl_mifd, cel_size, 'UniformOutput',false);
cgl_mifd = double(cell2mat(cgl_mifd'));

cgl_mafd = handles.Measurements.CellGolgi.AreaShape_MaxFeretDiameter;
cgl_mafd = cellfun(@(x,y) padarray(x,y-length(x),NaN,'post'), cgl_mafd, cel_size, 'UniformOutput',false);
cgl_mafd = double(cell2mat(cgl_mafd'));

cgl_area = double(cell2mat(cgl_area));
cgl_conv = cell2mat(cgl_conv);
cgl_avgI = cell2mat(cgl_avgI);
cgl_intI = cell2mat(cgl_intI);
cgl_eqdi = cell2mat(cgl_eqdi);
cgl_soli = cell2mat(cgl_soli);
cgl_mrad = cell2mat(cgl_mrad);
cgl_vrad = cell2mat(cgl_vrad);
cgl_miax = cell2mat(cgl_miax);
cgl_maax = cell2mat(cgl_maax);
cgl_orie = cell2mat(cgl_orie);
cgl_mdis = cell2mat(cgl_mdis);
cgl_vdis = cell2mat(cgl_vdis);

cel_rna_numb = double(cell2mat(handles.Measurements.Cells.Children_Maskedsirna_Count'));
cel_rna_maxI = cell2mat(handles.Measurements.Cells.Intensity_MaxIntensity_CropRNA');

cel_stat = [cel_stat nuc_labl nuc_area nuc_comp nuc_soli nuc_form ...
    nuc_perc nuc_locy nuc_locx nuc_miax nuc_maax nuc_mifd nuc_mafd cel_area cel_gol_numb ...
    cel_kif_intI cel_kif_avgI cel_kif_medI cel_kif_maxI cel_kif_minI cel_kif_loQI ...
    cel_kif_upQI cel_mafd cel_mifd cel_miax cel_maax cel_orie cel_mgl_avgA ...
    cel_rna_loQI cel_rna_minI cel_mgl_maxA mgl_intI mgl_maxI mgl_avgI cgl_area ...
    cgl_avgI cgl_intI cgl_comp cgl_soli cgl_form cgl_mrad cgl_vrad cgl_miax ...
    cgl_maax cgl_orie cgl_mifd cgl_mafd cgl_mdis cgl_vdis cel_rna_numb cel_rna_maxI ...
    mgl_100_area mgl_150_area mgl_225_area mgl_338_area mgl_per_area cgl_conv cgl_eqdi mgl_dis_area];
cel_legn = {'box' 'slid' 'slip' 'site' 'nuc_labl' 'nuc_area' 'nuc_comp' 'nuc_soli' 'nuc_form' ...
    'nuc_perc' 'nuc_locy' 'nuc_locx' 'nuc_miax' 'nuc_maax' 'nuc_mifd' 'nuc_mafd' ...
    'cel_area' 'cel_gol_numb' 'cel_kif_intI' 'cel_kif_avgI' 'cel_kif_medI' 'cel_kif_maxI' ...
    'cel_kif_minI' 'cel_kif_loQI' 'cel_kif_upQI' 'cel_mafd' 'cel_mifd' 'cel_miax' ...
    'cel_maax' 'cel_orie' 'cel_mgl_avgA' 'cel_rna_loQI' 'cel_rna_minI' 'cel_mgl_maxA' ...
    'mgl_intI' 'mgl_maxI' 'mgl_avgI' 'cgl_area' 'cgl_avgI' 'cgl_intI' ...
    'cgl_comp' 'cgl_soli' 'cgl_form' 'cgl_mrad' 'cgl_vrad' 'cgl_miax' 'cgl_maax' ...
    'cgl_orie' 'cgl_mifd' 'cgl_mafd' 'cgl_mdis' 'cgl_vdis' 'cel_rna_numb' 'cel_rna_maxI' ...
    'mgl_100_area' 'mgl_150_area' 'mgl_225_area' 'mgl_338_area' 'mgl_per_area' 'cgl_conv' 'cgl_eqdi' 'mgl_dis_area'}';