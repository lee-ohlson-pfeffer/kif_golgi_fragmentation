function [noc_mat_leg, noc_cel, noc_mat_exp, noc_exp, noc_cond] = postcp_noc_analysis(fileoutmat, ...
    noc_key, mean_kif_lo_thresh, mean_kif_hi_thresh, nuc_area_hi_thresh)

% post_noc_analysis take a *.mat file produced by cellprofiler
% (http://www.cellprofiler.org/) and extracts features, applies threshold
% to produce a figure of the Golgi fragmentation metric

% filematout is the *.mat file produced by cellprofiler

% noc_key is an mx7 matrix with m experiements descirbed by the 7 variables: 
% 'Experiment','Slide','Slip','Internal Exp Number','siRNA or Noc
% presence','Rescue Construct','Magnification' used to parse the different
% conditions.

% thresh inputs are mx1 arrays with experiment specific thresholds for
% cells that are non-rescued, too highly expressed, or large nuclei (likely
% a dividing cell.

% noc_mat_leg is a cell array describing the different metrics produced

% noc_cel is a mx2 cell array summarizing the results of m experiments with
% the legend for each experiment at {m,2}

% noc_mat_exp is a matrx with all results

% noc_exp is is a 1xn cell with the cummulated results from n experiements
% after thresholding each cell element is a 2x7 matrix representing the
% results for each condition within the experiment.

% noc_cond is a 2x7 cell with the results from noc_exp sorted into
% conditions. Each element contains an mx72 matrix with metrics for each of
% the m cells measured.

%% Initialization
[noc_mat_exp,noc_mat_leg] = postcp_noc(fileoutmat);


noc_cel = cell(max(noc_key(:,4)),2);

for for0 = 1:max(noc_key(:,4))
    [noc_key_exp,~] = find(noc_key(:,4)==for0);
    
    noc_mat_exp_sub = [];
    
    for for1 = 1:numel(noc_key_exp)
        
        [dif_cells,~] = find(noc_mat_exp(:,1)==noc_key(noc_key_exp(for1),1) ...
            & noc_mat_exp(:,2)==noc_key(noc_key_exp(for1),2) ...
            & noc_mat_exp(:,3)==noc_key(noc_key_exp(for1),3));
        
        noc_mat_exp_sub = cat(1,noc_mat_exp_sub,noc_mat_exp(dif_cells,:));
    end
    
    noc_cel{for0,1} = noc_mat_exp_sub;
    noc_cel{for0,2} = noc_key(noc_key_exp,:);
    
end

%% Thresholds

noc_exp = {};

for forX = 1:max(noc_key(:,4))
    exp_num_1 = forX;
    
    Adata = noc_cel{exp_num_1,1};
    Akey = noc_cel{exp_num_1,2};
    
    %GOL Presence
    Atemp = sum(Adata,2);
    Adata(isnan(Atemp),:) = [];
    
    %NUC Percentage
    Adata = thresholdstatistics(Adata,-10,20);
    
    %NUC Size
    Adata = thresholdstatistics(Adata,6,0.5e4);
    Adata = thresholdstatistics(Adata,-6,nuc_area_hi_thresh(forX));
    
    %NUC Solidity
    Adata = thresholdstatistics(Adata,8,0.82);
    
    %NUC FormFactor
    Adata = thresholdstatistics(Adata,9,0.34);
    
    %Rescue
    
    [Apos,~] = find(Akey(:,6)>0);
    Amat_pos = [];
    
    for for1 = 1:numel(Apos)
        
        [dif_cells,~] = find(Adata(:,1)==Akey(Apos(for1),1) ...
            & Adata(:,2)==Akey(Apos(for1),2) ...
            & Adata(:,3)==Akey(Apos(for1),3));
        
        Amat_pos = cat(1,Amat_pos,Adata(dif_cells,:));
    end
    
    [Aneg,~] = find(Akey(:,6)==0);
    Amat_neg = [];
    
    for for1 = 1:numel(Aneg)
        
        [dif_cells,~] = find(Adata(:,1)==Akey(Aneg(for1),1) ...
            & Adata(:,2)==Akey(Aneg(for1),2) ...
            & Adata(:,3)==Akey(Aneg(for1),3));
        
        Amat_neg = cat(1,Amat_neg,Adata(dif_cells,:));
    end
    
    if isempty(Amat_neg)
        
        temp_slope1_pos = Amat_pos(:,21)./Amat_pos(:,20);
        
        goodIndx1 = find((temp_slope1_pos)>0.6);
        Amat_pos = Amat_pos(goodIndx1,:); %#ok<*FNDSB>
        
        Amat_pos = thresholdstatistics(Amat_pos,20,mean_kif_lo_thresh(forX));
        Amat_pos = thresholdstatistics(Amat_pos,-20,mean_kif_hi_thresh(forX));
        
    else
        
        temp_slope1_pos = Amat_pos(:,21)./Amat_pos(:,20);
        
        goodIndx1 = find((temp_slope1_pos)>0.6);
        Amat_pos = Amat_pos(goodIndx1,:);
        
        Amat_pos = thresholdstatistics(Amat_pos,20,mean_kif_lo_thresh(forX));
        Amat_pos = thresholdstatistics(Amat_pos,-20,mean_kif_hi_thresh(forX));
        
    end
    
    Adata = cat(1,Amat_neg, Amat_pos);
    
    % Golgi weights
    
    % Column 70
    Amdis_weight = quantile(Adata(:,42),[0.2])./Adata(:,42); %#ok<*NBRAK>
    Adata = [Adata Amdis_weight];
    
    % Column 71
    Amdis_weight = quantile(Adata(:,49),[0.2])./Adata(:,49);
    Adata = [Adata Amdis_weight];
    
    % Column 72 normalize based on rescue intensity
    [A_exp,~] = find(Akey(:,5)==1&Akey(:,6)==1);
    
    noc_mat_exp = [];
    
    for for2 = 1:numel(A_exp)
        
        [dif_cells,~] = find(Adata(:,1)==Akey(A_exp(for2),1) ...
            & Adata(:,2)==Akey(A_exp(for2),2) ...
            & Adata(:,3)==Akey(A_exp(for2),3));
        
        noc_mat_exp = cat(1,noc_mat_exp,Adata(dif_cells,:));
    end
    
    Amdis_weight_mod = mean(noc_mat_exp(:,20));
    
    Amdis_weight = Adata(:,20)./Amdis_weight_mod;
    [A_exp,~] = find(Akey(:,6)<1);
    
    noc_mat_exp = [];
    for for2 = 1:numel(A_exp)
        [dif_cells,~] = find(Adata(:,1)==Akey(A_exp(for2),1) ...
            & Adata(:,2)==Akey(A_exp(for2),2) ...
            & Adata(:,3)==Akey(A_exp(for2),3));
        
        noc_mat_exp = cat(1,noc_mat_exp,dif_cells);
    end
    
    Amdis_weight(noc_mat_exp) = 1;
    Adata = [Adata Amdis_weight]; %#ok<*AGROW>
    
    
    %Split Conditions
    Acel = cell(max(Akey(:,5))+1,7);
    
    for for0 = 0:max(Akey(:,5))
        
        for for1 = 0:6
            
            [A_exp,~] = find(Akey(:,5)==for0&Akey(:,6)==for1);
            
            noc_mat_exp = [];
            
            if isempty(A_exp)
                
                noc_mat_exp = NaN(1,72);
                
            else
                
                for for2 = 1:numel(A_exp)
                    
                    [dif_cells,~] = find(Adata(:,1)==Akey(A_exp(for2),1) ...
                        & Adata(:,2)==Akey(A_exp(for2),2) ...
                        & Adata(:,3)==Akey(A_exp(for2),3));
                    
                    noc_mat_exp = cat(1,noc_mat_exp,Adata(dif_cells,:));
                end
                
            end
            
            
            Acel{for0+1,for1+1} = noc_mat_exp;
        end
    end
    
    noc_exp{exp_num_1} = Acel;
    
end

%% Graphing

noc_cond = cell(2,7);

for for1 = 1:numel(noc_exp)
    temp_exp = noc_exp{for1};
    for for2 = 1:size(temp_exp,1)
        for for3 = 1:size(temp_exp,2)
            temp_cond = temp_exp{for2,for3};
            noc_cond{for2,for3} = cat(1,noc_cond{for2,for3},temp_cond);
        end
    end 
end

num1 = 51;
div1 = 72;
figure
temp_norm = cellfun(@(x) x(:,num1)./x(:,div1),noc_cond,'UniformOutput',0);
a = cellfun(@(x) nanmean(x),temp_norm);
b = cellfun(@(x) nanstd(x),temp_norm);
c = cellfun(@(x) sum(~isnan(x)),temp_norm);
d = c.^0.5;
e = b./d;
a1 = [a(1,:) a(2,:)];
e1 = [e(1,:) e(2,:)];

barweb(a1,e1)

