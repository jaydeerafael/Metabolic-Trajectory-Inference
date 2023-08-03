% 1. Load single-cell data.

filename = '/Users/jaydee/Desktop/Rannís_proj/GSE123813_bcc_scRNA_counts.txt'; % Replace filepath

% Read data without considering the first row as variable names
data_bcc = readtable(filename, 'Delimiter', '\t', 'ReadVariableNames', false);

% Rename the first column into Gene_IDs and the rest as "sample 1", and etc.
data_bcc.Properties.VariableNames{1} = 'Gene_IDs';
numOfColumns = size(data_bcc,2);
sampleNames = arrayfun(@(x) sprintf('sample %d', x), 1:numOfColumns-1, 'UniformOutput', false);
data_bcc.Properties.VariableNames(2:end) = sampleNames;
%% 
% 
% 2. Separate data into two dataframes.

gene_ids = data_bcc.Gene_IDs;
matrix_data = table2array(data_bcc(:, 2:end));
%% 
% 
% 2a. Make a smaller dataset for testing and prodding.

% Selecting all the rows/genes but only 10 or x many cells/samples you want
smaller_data = data_bcc(:, 1:10);

% separate smaller_data into two dataframes 
gene_ids = smaller_data.Gene_IDs;
smaller_matrix_data = table2array(smaller_data(:, 2:end));
%% 
% 
% 3. Create a column name dataframe.

% Number of samples
n_samples = size(data_bcc, 2);

% Create cell array with specified names
colnames_data = cell(1, n_samples);

% for smaller_data
n_samples = size(smaller_data, 2);
smaller_colnames_data = cell(1, n_samples);
%% 
% 
% 4. initiate Cobra Toolbox as a solver:

addpath(genpath('/Applications/CPLEX_Studio129/cplex/matlab/x86-64_osx/')) % Replace filepath
addpath(genpath('/Users/jaydee/cobratoolbox'))  % Replace filepath
initCobraToolbox(false)
%% 
% 
% 5. Load Dictionary and Recon3D Model.

load('dico_rFASTCORMICS.mat')  % dictionary to map the rownname identifier to the genes in the model
load('/Users/jaydee/cobratoolbo/test/verifiedTests/visualization/Recon3D_301/Recon3DModel_301.mat') % Replace filepath
model = Recon3DModel
%% 
% 
% 6. Perform a Consistency Check (FastCC) on Recon3D.

A_fast = fastcc_4_rfastcormics(Recon3DModel, 1e-4,0) % Create consistent model by running FASTCC (Vlassis et al., 2014)
Cmodel_fast = removeRxns(Recon3DModel, Recon3DModel.rxns(setdiff(1:numel(Recon3DModel.rxns),A_fast))) % Removing blocked reactions
%% 
% 
% 7. Check Objective Function (biomass_maintenance) of Recon3D and Cmodel_fast.

objectiveRecon3D = checkObjective(model)
objectiveCmodel = checkObjective(Cmodel_fast) 
biomass_rxn = checkObjective(Cmodel_fast) % Rename the biomass reaction innto biomass_rxn as a separate variable in workspace
printRxnFormula(model, objectiveCmodel) % Too basically "confirm" that the model does infact, have a biomass reaction
%% 
% 
% 8. Optional Inputs.

already_mapped_tag = 0;
epsilon = 1e-4; % Avoid small number errors
consensus_proportion = 0.9; % Gene has to be expressed in 90% of the cases in order to be included. Only relevant if you want to create one generic model from different samples

optional_settings.func = {'biomass_maintenance','DM_atp_c_'};
not_medium_constrained = '0'
optional_settings.not_medium_constrained = not_medium_constrained;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Other optional inputs if more is known of the medium
load medium_example % need to define medium for the cells used here

unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};
unpenalized = Cmodel_fast.rxns(ismember(Cmodel_fast.subSystems,unpenalizedSystems));

optional_settings.unpenalized = unpenalized;
optional_settings.func = {'biomass_reaction','DM_atp_c_'}; % Forced additional reactions into the  model
not_medium_constrained = 'EX_tag_hs(e)';
optional_settings.not_medium_constrained = not_medium_constrained;
optional_settings.medium = medium_example; % Remove field if no constraint is provided
%% 

% 9. Create Context Specific Models.

for i = 1:numel(colnames_data) % For each sample
    [model_out{i}, A_keep{i}] = fastcormics_RNAseq(Cmodel_fast, matrix_data(:,i), gene_ids, dico, ...
        biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
end
%% 

% 10. Perform FBA and pFBA on the models

% FBA
% Number of models in model_out
numModels = numel(model_out);

% Initialize an empty cell array to store FBA scores
FBAscores = cell(numModels, 1);

% Loop through each model and perform FBA
for i = 1:numModels
    result = optimizeCbModel(model_out{i});
    FBAscores{i} = result;  % Assuming result is a structure. If it's a simple value, you can use an array instead of a cell array.
end

% Convert the cell array to a table
FBAscoreTable = cell2table(FBAscores, 'VariableNames', {'FBA_Result'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pFBA
% Number of models in model_out
numModels = numel(model_out);

% Initialize an empty cell array to store pFBA results
pFBAscores = cell(numModels, 1);

% Loop through each model and perform pFBA
for i = 1:numModels
    result = optimizeCbModel(model_out{i}, 'parsimoniousFBA');
    pFBAscores{i} = result;  % Assuming result is a structure. If it's a simple value, you can use an array instead of a cell array.
end

% Convert the cell array to a table
pFBAscoreTable = cell2table(pFBAscores, 'VariableNames', {'pFBA_Result'});

