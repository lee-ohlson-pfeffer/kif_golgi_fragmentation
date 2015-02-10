function [rna_mat_leg, rna_cel, rna_mat_exp, rna_exp, rna_cond] = postcp_sirna_analysis(fileoutmat, ...
    rna_key, max_rna_lo_thresh, mean_kif_lo_thresh, mean_kif_hi_thresh)

% post_sirna_analysis take a *.mat file produced by cellprofiler
% (http://www.cellprofiler.org/) and extracts features, applies threshold
% to produce a figure of the Golgi fragmentation metric

% filematout is the *.mat file produced by cellprofiler

% rna_key is an mx7 matrix with m experiements descirbed by the 7 variables: 
% 'Experiment','Slide','Slip','Internal Exp Number','siRNA or Noc
% presence','Rescue Construct','Magnification' used to parse the different
% conditions.

% thresh inputs are mx1 arrays with experiment specific thresholds for
% cells that are non-rescued, too highly expressed, or large nuclei (likely
% a dividing cell.

% rna_mat_leg is a cell array describing the different metrics produced

% rna_cel is a mx2 cell array summarizing the results of m experiments with
% the legend for each experiment at {m,2}

% rna_mat_exp is a matrx with all results

% rna_exp is is a 1xn cell with the cummulated results from n experiements
% after thresholding each cell element is a 2x7 matrix representing the
% results for each condition within the experiment.

% rna_cond is a 2x7 cell with the results from rna_exp sorted into
% conditions. Each element contains an mx72 matrix with metrics for each of
% the m cells measured.

%% Initialization
[rna_mat_exp,rna_mat_leg] = postcp_sirna(fileoutmat);


rna_cel = cell(max(rna_key(:,4)),2);

for for0 = 1:max(rna_key(:,4))
    [rna_key_exp,~] = find(rna_key(:,4)==for0);
    
    rna_mat_exp_sub = [];
    
    for for1 = 1:numel(rna_key_exp)
        
        [dif_cells,~] = find(rna_mat_exp(:,1)==rna_key(rna_key_exp(for1),1) ...
            & rna_mat_exp(:,2)==rna_key(rna_key_exp(for1),2) ...
            & rna_mat_exp(:,3)==rna_key(rna_key_exp(for1),3));
        
        rna_mat_exp_sub = cat(1,rna_mat_exp_sub,rna_mat_exp(dif_cells,:));
    end
    
    rna_cel{for0,1} = rna_mat_exp_sub;
    rna_cel{for0,2} = rna_key(rna_key_exp,:);
    
end

%% Thresholds
rna_exp = {};

for forX = 1:max(rna_key(:,4))
    exp_num_1 = forX;
    
    Adata = rna_cel{exp_num_1,1};
    Akey = rna_cel{exp_num_1,2};
    
    %GOL Presence
    Atemp = sum(Adata,2);
    Adata(isnan(Atemp),:) = [];
    
    %NUC Percentage
    Adata = thresholdstatistics(Adata,-10,20);
    
    %NUC Solidity
    Adata = thresholdstatistics(Adata,8,0.84);
    
    %NUC FormFactor
    Adata = thresholdstatistics(Adata,9,0.35);
    
    %GOL Area
    Adata = thresholdstatistics(Adata,38,500);
    
    %RNA
    [Apos,~] = find(Akey(:,5)==1);
    Amat_pos = [];
    
    for for1 = 1:numel(Apos)
        
        [dif_cells,~] = find(Adata(:,1)==Akey(Apos(for1),1) ...
            & Adata(:,2)==Akey(Apos(for1),2) ...
            & Adata(:,3)==Akey(Apos(for1),3));
        
        Amat_pos = cat(1,Amat_pos,Adata(dif_cells,:));
    end
    
    [Aneg,~] = find(Akey(:,5)~=1);
    Amat_neg = [];
    
    for for1 = 1:numel(Aneg)
        
        [dif_cells,~] = find(Adata(:,1)==Akey(Aneg(for1),1) ...
            & Adata(:,2)==Akey(Aneg(for1),2) ...
            & Adata(:,3)==Akey(Aneg(for1),3));
        
        Amat_neg = cat(1,Amat_neg,Adata(dif_cells,:));
    end
    
    Amat_pos = thresholdstatistics(Amat_pos,54,max_rna_lo_thresh(forX));
    
    Adata = cat(1,Amat_neg, Amat_pos);
    
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
    
    
    
    temp_slope1_pos = Amat_pos(:,21)./Amat_pos(:,20);
    goodIndx1 = find((temp_slope1_pos)>0.6);
    Amat_pos = Amat_pos(goodIndx1,:); %#ok<*FNDSB>
    
    
    Amat_pos = thresholdstatistics(Amat_pos,20,mean_kif_lo_thresh(forX));
    Amat_pos = thresholdstatistics(Amat_pos,-20,mean_kif_hi_thresh(forX));
    
    Adata = cat(1,Amat_neg, Amat_pos);
    
    % Golgi weights
    
    % Col 74
    Amdis_weight = (Adata(:,42)+Adata(:,59)+quantile(Adata(:,44),[0.2])./Adata(:,44))./3; %#ok<*NBRAK>
    Adata = [Adata Amdis_weight]; %#ok<*AGROW>
    
    % Col 75
    Amdis_weight = quantile(Adata(:,44),[0.2])./Adata(:,44);
    Adata = [Adata Amdis_weight];
    
    % Col 76
    Amdis_weight = Adata(:,57);
    Adata = [Adata Amdis_weight];
    
    % Col 77 normalize based on wild-type rescue intensity
    [A_exp,~] = find(Akey(:,5)==1&Akey(:,6)==1);
    
    rna_mat_exp = [];
    
    for for2 = 1:numel(A_exp)
        
        [dif_cells,~] = find(Adata(:,1)==Akey(A_exp(for2),1) ...
            & Adata(:,2)==Akey(A_exp(for2),2) ...
            & Adata(:,3)==Akey(A_exp(for2),3));
        
        rna_mat_exp = cat(1,rna_mat_exp,Adata(dif_cells,:));
    end
    
    Amdis_weight_mod = mean(rna_mat_exp(:,20));
    
    Amdis_weight = Adata(:,20)./Amdis_weight_mod;
    [A_exp,~] = find(Akey(:,6)<1);
    
    rna_mat_exp = [];
    for for2 = 1:numel(A_exp)
        [dif_cells,~] = find(Adata(:,1)==Akey(A_exp(for2),1) ...
            & Adata(:,2)==Akey(A_exp(for2),2) ...
            & Adata(:,3)==Akey(A_exp(for2),3));
        
        rna_mat_exp = cat(1,rna_mat_exp,dif_cells);
    end
    
    Amdis_weight(rna_mat_exp) = 1;
    Adata = [Adata Amdis_weight];
    
    % Col 78
    Amdis_weight = Adata(:,end-1)./Adata(:,end);
    Adata = [Adata Amdis_weight];
    
    % Col 79
    Amdis_weight_mod = 0.0032706;
    Amdis_weight = Adata(:,20)./Amdis_weight_mod;
    Amdis_weight(rna_mat_exp) = 1;
    Adata = [Adata Amdis_weight];
    
    % Col 80 = (Col 76 = Col 57)./Col 79 - Percent Golgi in large structures
    % divided by the normalized intensity
    Amdis_weight = Adata(:,end-3)./Adata(:,end);
    Adata = [Adata Amdis_weight];
    
    % Col 81
    Amdis_weight = Adata(:,42)./Adata(:,79);
    Adata = [Adata Amdis_weight];
    
    % Col 82
    Amdis_weight = Adata(:,75)./Adata(:,79);
    Adata = [Adata Amdis_weight];
    
    
    %Split Conditions
    Acel = cell(max(Akey(:,5)),max(Akey(:,6)));
    
    for for0 = 0:max(Akey(:,5))
        
        for for1 = 0:max(Akey(:,6))
            
            [A_exp,~] = find(Akey(:,5)==for0&Akey(:,6)==for1);
            
            rna_mat_exp = [];
            
            if isempty(A_exp)
                
                rna_mat_exp = NaN(1,82);
                
            else
                
                for for2 = 1:numel(A_exp)
                    
                    [dif_cells,~] = find(Adata(:,1)==Akey(A_exp(for2),1) ...
                        & Adata(:,2)==Akey(A_exp(for2),2) ...
                        & Adata(:,3)==Akey(A_exp(for2),3));
                    
                    rna_mat_exp = cat(1,rna_mat_exp,Adata(dif_cells,:));
                end
                
            end
            
            
            Acel{for0+1,for1+1} = rna_mat_exp;
            
        end
    end
    
    rna_exp{exp_num_1} = Acel;
    
end

%% Graphing

rna_cond = cell(2,7);

for for1 = 1:numel(rna_exp)
    temp_exp = rna_exp{for1};
    for for2 = 1:size(temp_exp,1)
        for for3 = 1:size(temp_exp,2)
            temp_cond = temp_exp{for2,for3};
            rna_cond{for2,for3} = cat(1,rna_cond{for2,for3},temp_cond);
        end
    end
end

for1 = 80;
figure
Col1 = for1;
a = cellfun(@(x) nanmean(x(:,Col1)),rna_cond);
b = cellfun(@(x) nanstd(x(:,Col1)),rna_cond);
c = cellfun(@(x) sum(~isnan(x(:,Col1))),rna_cond);
d = c.^0.5;
e = b./d;
a1 = [a(1,:) a(2,:)];
e1 = [e(1,:) e(2,:)];

barweb(a1,e1)
