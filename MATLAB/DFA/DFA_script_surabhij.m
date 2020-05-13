%Step 1: make sure COBRA is installed and it works 
    initCobraToolbox;
    changeCobraSolver('gurobi');
%Step 2: DFA input - Metabolic model - make it into biomass objective
    load('CoreModel.mat');
    metabolicmodel = changeObjective(core_genecomb, 'biomass_NCI60');
%Step 3: DFA input - metabolomic - create struct with position and data
    %3a: read in excel sheet and turn into table
   core_metabolomics = readtable('tutorial.xlsx', ...
    'Sheet', 'CORE', ...
    'Format','auto', ...
    'ReadRowNames', true);
    %3b: save columns, positions and data
    columns = core_metabolomics.Properties.VariableNames;
    positions = table2array(core_metabolomics(:, 'positionInModel'));
    a498_data = table2array(core_metabolomics(:, contains(columns, 'A498')));
    loximvi_data = table2array(core_metabolomics(:, contains(columns, 'LOXIMVI')));
    %3c: create struct
    a498_metabolomics.positions    = positions; 
    loximvi_metabolomics.positions = positions;
    a498_metabolomics.data         = a498_data;      
    loximvi_metabolomics.data      = loximvi_data;
%Step 4: set parameteres
    params.kappa  = 1E-1;
    params.kappa2 = 1E-6;
    params.norm   = 'None';
%Step 5: Run DFA
    [a498_model, a498_soln] = DFA(metabolicmodel, a498_metabolomics,params);
    [loximvi_model, loximvi_soln] = DFA(metabolicmodel, loximvi_metabolomics,params);
    
%Step 6: Create Heatmap of reactions from DFA model and biomass solution.x 
    figure (1);
    val = cell2table([a498_model.rxns,num2cell(a498_soln.x)]);
    heatmap(val);
    