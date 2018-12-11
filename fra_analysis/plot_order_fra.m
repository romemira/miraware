%% Define params
fra_psth_path = '/media/miranda/09F6A5BC59F7E39C/metadata/fra_psth/Ronnie_P05_quning_normal'; %points to fra_psth file
save_dir = '/media/miranda/09F6A5BC59F7E39C/FRA_plots/Ronnie/P05_normal/'; %Points to where we wanna save the plots
p_val_thr = 0.0001; %The pvalue threshold that we sue to select FRAs

%% Load and sort
load(fra_psth_path);
temp_p = fra_psth(1).pval; %Create temp var to work on 
temp_p_sort = sortrows(temp_p,1,'ascend'); %Sort based on the pval for freqeuncy from anova test 
cluster_id_sort = temp_p_sort(:,4); %Get the sorted ids of all clusters
p_val_freq_sort = temp_p_sort(:,1); %Get the sorted freq p vals for all clsuters
clusters_id_final = cluster_id_sort(p_val_freq_sort<p_val_thr); %Select only clusters that pass the freq p val threshold
%% Run function
plot_fra_pixels(fra_psth,clusters_id_final,save_dir);