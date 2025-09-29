function cationCalc4epma_eng_250925()
% tested with matlab 2024b
% Kazuki Matsuyama (29/08/2025 @ University of Montpellier, France)
% Yumiko Harigane (@ GSJ, AIST, Japan)
% Yoshihiro Nakamura (@ GSJ, AIST, Japan)
%
% Syntax:
% function cationCalc4epma_250827()
%
% Inputs:
% ------
% - None
%
% Outputs:
% ------
% - None
%

%% 0.Clear workspace
% Deleting variables (cleaning workspace)
clear all

% Closing windows
close all

% Measuring time
tic

%% 1.Import Excel files
% Load the Excel file (raw data) for performing the cation calculation and the Excel file (Cation_module) containing the modules used in the calculation into the workspace.
% The raw data is assumed to be in .csv or .xlsx format.
disp('------------------------------------------------------------------');
disp('=== Step 1: Importing Excel files ===');

% select a .csv or .xlsx file of raw data using list dialog
[excel_dataset, ~] = uigetfile( ...
    {'*.xlsx;*.csv','Excel or CSV files (*.xlsx, *.csv)'}, ...
    'Select .csv or .xlsx file', ...
    'MultiSelect', 'off');

% If an excel file is NOT selected, it will stop with an error
if isequal(excel_dataset, 0)
    error('File selection was canceled.');
end

% If the file name is already determined, user can specify the file name directly
% excel_dataset = 'EPMA_test.xlsx';

% Display the file name of the imported raw data
disp([excel_dataset ' has been imported as raw data.' ]);

% Acquire the original file name (excluding the extension)
[~, filename, ~] = fileparts(excel_dataset);

% Load EPMA dataset into the workspace
% The leftmost sheet is automatically selected
imported_data = readtable(excel_dataset);

% Load the cation calculation module (.xlsx file)
disp('Reading cation calculation dataset...');

% Specify the name of the module
excel_cation = 'Cation_moduli.xlsx';
% Load 'Molar weight' sheet
moduli_molarWeight = readcell(excel_cation, 'Sheet', 'Molar weight');

% Load 'Stoichiometry' sheet
moduli_Stoichiometry = readtable(excel_cation, 'Sheet', 'Stoichiometry');

% Report loading completion
disp('=== Cation calculation dataset has been imported ===');

%% 2.Select mineral phases
% Select target minerals from the list in "Cation_moduli".
% Create a list of mineral names to use for mineral phase identification (section 9).
% This script performs uniform processing based on all contents in this list.
% If users want to add their own target minerals, 
% add the information to the 'Stoichiometry' sheet in Cation_moduli.xlsx.
% The conditions for mineral phase identification are set in the section 9.
disp('------------------------------------------------------------------');
disp('=== Step 2: Defining mineral phase list ===');

% Select mineral phases using list dialog
[minerals_id, ok] = listdlg( ...
    'PromptString', 'Select mineral phases for cation calculation:', ...
    'ListString', moduli_Stoichiometry.Properties.VariableNames(2:end), ...
    'SelectionMode', 'multiple', ...
    'ListSize', [300 300], ...
    'InitialValue',1:numel(moduli_Stoichiometry.Properties.VariableNames)-1);

if ok == 0
    % If any phases are NOT selected, it will stop with an error
    error('No mineral phases selected. Operation canceled.');
end

% Make a list of mineral names
minerals = moduli_Stoichiometry.Properties.VariableNames(minerals_id+1);

% If users want to process uniformly, they can specify it here:
% In that case, comment out the list dialog section above
%minerals = ["Olivine", "Opx", "Cpx", "Spinel", "Plagioclase", "Ilmenite", "Magnetite", "Quartz", "Amphibole", "Chlorite", "Epidote", "Mica", "Garnet"];

% Display the list of mineral names
disp(minerals);

%% 3.Pre-process the imported dataset
% Before starting the cation calculation, preprocess the raw data.
disp('------------------------------------------------------------------');
disp('=== Step 3: Pre-processing Excel files ===');

% Exclude data with missing comments from imported EPMA data
% Delete rows where the "Comment" column is NaN
% NB: the existence of the "Comment" column is implicitly assumed
if any(strcmp('Comment', imported_data.Properties.VariableNames))
    % Create a logical index: find rows with NaN
    nan_rows = cellfun(@(x) isempty(x) || strcmpi(x, 'NaN'), imported_data.Comment);

    % Display the number of rows with NaN
    disp(['Removing ', num2str(sum(nan_rows)), ' rows with missing Comment.']);

    % Delete all rows with NaN
    imported_data(nan_rows, :) = [];

    % If all the "Comment" columns in imported_data are empty,
    % it will stop with an error
    if isempty(imported_data.Comment)
        error('All "Comment" rows were removed during preprocessing. Check the original file.');
    end

else
    % If "Comment" column does not exist, it will stop with an error
    error('"Comment" column not found in the imported table.');
end

% Check the size of imported_data
% size_imported_data = [height, width]
size_imported_data = size(imported_data);
height_imported_data = size_imported_data(1);

% Report pre-processing completion
disp('=== Pre-processing has been finished ===');

%% 4.Select data to calculate cations
% From the raw data loaded, 
% select the data to which users want to apply the cation calculation.
% The selected data is assumed to be mass% or K-ratio value data.
% This script assumes that the selected data for each analysis point
% (each row) totals approx. 100 mass%.
% If the same element (or oxide) is selected,
% it will stop with an error along with a warning message.
% The data selected here will be output as "filename_selected.xlsx."
disp('------------------------------------------------------------------');
disp('=== Step 4: Selecting data columns for cation calculation ===');

% Create an empty table to store the selected data
selected_data = table();

% Insert the "Comment" column of the raw data into the leftmost column of selected_data
if any(strcmp('Comment', imported_data.Properties.VariableNames))
    selected_data.Comment = imported_data.Comment;
else
    % If "Comment" column does not exist, it will stop with an error
    error('"Comment" column not found in Table_data.');
end

% All variable names in imported_data are displayed in a list dialog (except for "Comment")
variable_list = imported_data.Properties.VariableNames;
% Exclude "Comment" column
variable_list(strcmp(variable_list, 'Comment')) = []; 

% Select data to calculate cations using list dialog
[selected_vars, ok] = listdlg( ...
    'PromptString', 'Select data for cation calculation:', ...
    'ListString', variable_list, ...
    'SelectionMode', 'multiple', ...
    'ListSize', [300 300]);

if ok == 0
    % If any data is not selected, it will stop with an error
    error('No variables selected. Operation canceled.');
end

% Store the selected data in selected_data along with the original variable names
selected_variable_names = variable_list(selected_vars);
for i = 1:length(selected_variable_names)
    selected_data.(selected_variable_names{i}) = imported_data.(selected_variable_names{i});
end

% Check if the same element is selected (except for "Comment")
var_names = setdiff(selected_data.Properties.VariableNames, {'Comment'});

% Extract the beginning of the variable name (before the first "_")
% NB: assuming it is the names of elements or oxides
% If they contain no "_", proceed as is
oxides_names = cellfun(@(x) strtok(x, '_'), var_names, 'UniformOutput', false);

% Check for duplicates
[unique_names, ~, idx] = unique(oxides_names);
duplicate_elements = unique_names(accumarray(idx,1) >= 2);

% If there are duplicates, it will stop with an error
if ~isempty(duplicate_elements)
    error(['The same element is selected: ', strjoin(duplicate_elements, ', ')]);
end

% Export the table to an external file 
% with the name "original file name_selected.xlsx"
output_filename = strcat(filename, '_Selected.xlsx');
% Set the sheet name to "Selected raw data"
try
    writetable(selected_data, output_filename, 'Sheet','Selected raw data');
    disp(['Selected raw data written to: ', output_filename]);
catch ME
    error(['Failed to write output file: ', ME.message]);
end

% Report selecting data completion
disp('=== Selecting data columns has been finished ===');

%% 5.Apply Detection Limit (DL) filtering
% If the raw data contains Detection Limit (DL),
% DL filtering is performed on the selected data set.
% Specifically, if any data falls below the DL value, 
% that value is replaced with NaN.
% NB: In subsequent calculations, NaN is treated as 0:
% but NaN replacement is used to clarify whether filtering has occurred.
% DL and raw data sre assumed to be expressed in ppm and %, respectively:
% therefore the DL value is multiplied by 10^-4.
% Filtering is performed using the local function (1) apply_DL_mask.
disp('------------------------------------------------------------------');
disp('=== Step 5: Applying Detection Limit (DL) filtering ===');

% Extract columns containing "D_L" from imported_data
DL_columns = contains(imported_data.Properties.VariableNames, 'D_L');

% If the DL column does NOT exist, skip this process.
if ~any(DL_columns)
    disp('No detection limit (D_L) columns found in imported_data. Skipping DL filtering.');
    % If "D_L" column exists, perform filtering
else
    % Extract "D_L" column
    DLs = imported_data(:, DL_columns);

    % Define a function to extract the oxide name (prefix)
    extract_prefix = @(names) cellfun(@(x) strtok(x, '_'), names, 'UniformOutput', false);

    % For the variables in the "D_L" column, extract the element or oxide names
    % Extract the variable names in the "D_L" column
    DL_names = DLs.Properties.VariableNames;
    % Create a cell array to store the new variable names
    new_DL_names = extract_prefix(DL_names);

    % Extract the string before the first "_" from the variable names in selected_data
    oxides_names = extract_prefix(selected_variable_names);

    % Mltiply DLs by 10^-4 (to convert from ppm to %)
    DLs_pct = DLs{:,:} * 1e-4;

    % Check oxides_names from the first one, 
    % and if the string before the first "_" matches the string in the D_L column, 
    % perform the DL process.
    for i = 1:length(oxides_names)
        % Search for the name of an oxide that matches 
        % exactly from the first line of moduli_molarWeight
        id1 = strcmp(oxides_names{i}, moduli_molarWeight(1,:));

        if ~any(id1)
            % If an oxide is not found, it will stop with an error
            error(['Oxide name not found in moduli_molarWeight: ', oxides_names{i}]);
        end

        % Extract the element name corresponding to the oxide name
        % (e.g., SiO2 => Si)
        element = moduli_molarWeight{5,id1};

        % Search for the position of the DLs that perfectly matches the extracted element name
        id2 = find(strcmp(element, new_DL_names));

        if isempty(id2)
            error(['DL does not exist within imported dataset: ', oxides_names{i}]);
        elseif numel(id2) > 1
            warning(['Multiple DLs matched: ', oxides_names{i}, ' → Use the first match.']);
            id2 = id2(1);
        end

        % Extract the corresponding DL value
        DL = DLs_pct(:,id2);

        % If the Mass% data is less than the DL value, replace the value with NaN
        % use the local function (1) apply_DL_mask
        selected_data_i = contains(selected_data.Properties.VariableNames, oxides_names{i});
        if any(selected_data_i) && ~isempty(DL)
            selected_data{:, selected_data_i} = apply_DL_mask( ...
                selected_data{:, selected_data_i}, DL);
        else
            % If there is no corresponding raw data or DL, 
            % it will stop with an error
            error('No columns found in selected_data or DL');
        end

    end
end

% Report DL filtering completion
disp('=== DL filtering has been finished ===');

%% 6.Setup variable names
% Rename the selected_data variables.
% Specifically, if the variable name in the original raw data table is not
% just the oxide name (e.g., SiO2_Mass%), delete the part after the "_".
% Exclude the "Comment" columns from this process.
% If there are duplicate variable names, the process will stop with an error.
% After renaming the variables,
% selected_data will be output as an Excel file at the end of the section.
disp('------------------------------------------------------------------');
disp('=== Step 6: Renaming variables by removing suffixes after "_" ===');

% Exclude special columns (those that do not change)
var_names = selected_data.Properties.VariableNames;
exclude_names = {'Comment'};

% Create a cell array to store the new variable names
new_var_names = var_names;
for i = 1:numel(var_names)
    varname = var_names{i};

    % Process variable names that are not in the exclusion list and contain "_"
    if ~ismember(varname, exclude_names) && contains(varname, '_')
        % Extract the part before the first "_"
        new_var_names{i} = strtok(varname, '_');
    end
end

% Check for duplicate variable names:
% If there are duplicates, it will stop with an error
if numel(unique(new_var_names)) ~= numel(new_var_names)
    error('Renaming would result in duplicate variable names. Please revise variable selection.');
end

% Apply new variable names to selected_data
selected_data.Properties.VariableNames = new_var_names;

% Export the table to an external file 
% with the name "original file name_selected.xlsx"
% The sheet name is "Proc DL"
try
    writetable(selected_data, output_filename, 'Sheet','Proc DL');
    disp(['Selected data written to: ', output_filename]);
catch ME
    error(['Failed to write output file (Proc DL sheet): ', ME.message]);
end

% Report exporting a file
disp(['Selected data written to: ', output_filename]);

% Report Renaming variables completion
disp('=== Renaming variables has been finished ===');

%% 7.Process elements NOT present in the measured data
% For elements NOT included in the measurement data,
% add a new column and store NaN.
% Possible elements are obtained from moduli_molarWeight.
disp('------------------------------------------------------------------');
disp('=== Step 7: Processing elements not present in the measurement data ===');

% Obtain possible elements through EPMA analysis
elements_candidate = string(moduli_molarWeight(1, 2:end));

% Check the size of selected_data
% size_selected_data = [height, width]
size_selected_data = size(selected_data);
height_selected_data = size_selected_data(1);
nans_selected_data = NaN(height_selected_data,1);

% For elements_candidate that are NOT included in selected_data,
% add variable names and assign NaN to them all
for i = 1:length(elements_candidate)
    if ~ismember(elements_candidate(i), selected_data.Properties.VariableNames(:))
        selected_data.(elements_candidate(i)) = nans_selected_data;
        disp(['Column "'  elements_candidate{i}  '" has been added into selected_data with NaN values']);
    end
end

% Report processing completion
disp('===  Processing elements not present has been finished ===');

%% 8.Calculate Total mass%
% Sum the raw data for each analysis point (each row),
% and calculate the total value (Mass%).
% The calculated total value is added to the rightmost column of selected_data.
% NB: NaN is handled as 0 during the calculation process.
% NB: If this process is not performed, any calculation that contains 
% even one NaN will result in all results being NaN.
disp('------------------------------------------------------------------');
disp('=== Step 8: Making "Total" column ===');

% Acquire column names except for "Comment" column
numeric_var_names = setdiff(selected_data.Properties.VariableNames, {'Comment'});

% Check the data type of each column and select only numeric columns
valid_numeric_vars = {};
for i = 1:numel(numeric_var_names)
    col = selected_data.(numeric_var_names{i});
    if isnumeric(col) || islogical(col)
        valid_numeric_vars{end+1} = numeric_var_names{i}; %#ok<AGROW>
    end
end

% Sum rows by row while handling NaN as 0
if ~isempty(valid_numeric_vars)
    numeric_data = selected_data{:, valid_numeric_vars};

    % If the Total column is selected, display a warning message.
    if ismember('Total', selected_data.Properties.VariableNames)
        warning('"Total" column already exists and will be overwritten.');
    end

    % For each row, sum all the data
    selected_data.Total = sum(numeric_data, 2,'omitnan');
else
    % If a numeric column is not found and the Total value cannot be calculated,
    % it will stop with an error
    error('No numeric columns found (excluding "Comment"). "Total" column not created.');
end

% Report Making "Total" column completion
disp('===  Making "Total" column has been finished ===');

%% 9.Identify the mineral phases
% Automatically identify mineral phases from measurement data 
% by setting empirical thresholds.
% Classify representative igneous and metamorphic rock minerals.
% NB: Candidates of mineral phases are stored in moduli_Stoichiometry.
% The identification results are stored in "Phase" column in selected_data.
%
% This script uses an if statement 
% to make it easier for individual users to add their own mineral phases.
% NB: This classifications  are based solely on empirical rules.
% NB: Users should verify that the results are consistent with their analysis.
%
%
% About automatic identification:
% Its accuracy depends on the order in which the identification is made,
% so the it will be conducted based on the following criteria:
%
% 1) Identify mineral phases with distinctive compositions
% (i.e. minerals are unlikely to overlap with silicates).
% *The order of judgement is:
% phosphates => metal oxides => single components => low-Si oxides
%
% 2) Identify silicates in order of alkalinity, calcium, and water content.
% *Feldspar => Mica => Epidote => Chlorite
%
% 3) Identify Olivine/Pyroxenes
% Olivine and pyroxenes have similar compositions
% *Identify Olivine first (Si- and Ca-poor) to prevent misclassification
% **Then Cpx => Opx using the amount of CaO
%
% 4) Identify Garnet 
% (poor in alkali, and unique composition rich in Fe+Mn+Mg+Ca).
%
% 5) Identify Amphibole 
% (rich in alkali, Ca, OR Fe+Mg).
%

disp('------------------------------------------------------------------');
disp('=== Step 9: Mineral phases determination ===');

% NB：All thresholds are empirical:
% therefpre it is recommended to modify them according to the situation.
warning('Since the classification criteria are empirical, it is advisable to adjust them to fit the specific context.')

% Repeat the identification process starting from the top of Table_data
for i = 1:height_imported_data
    % "Excluded"
    % 1) total value <= 80 OR
    % 2) 102 <= total value
    if selected_data.Total(i) <= 80 || selected_data.Total(i) >= 102 % 80<=Total<=102
        % Phase naming "Excluded"
        selected_data.Phase(i) = "Excluded";

        % from here on, identification mineral phases
        % "Apatite"
        % 1) SiO2 <= 5 AND
        % 2) 50 <= CaO <= 60 AND
        % 3) 35 <= P2O5 <= 45
    elseif ismember('Apatite',minerals) && ...
            (selected_data.SiO2(i) <= 1 || isnan(selected_data.SiO2(i))) && ... % SiO2<=5
            (selected_data.CaO(i) >= 50 && selected_data.SiO2(i) <= 60) && ... % 50<=CaO<=60
            (selected_data.P2O5(i) >= 35 && selected_data.P2O5(i) <= 45) % 35<=P2O5<=45
        % Phase naming "Apatite"
        selected_data.Phase(i) = "Apatite";

        % "Ilmenite"
        % 1) SiO2 <= 5 AND
        % 2) 30 <= TiO2 AND
        % 3) 20 <= FeO
    elseif ismember('Ilmenite',minerals) && ...
            (selected_data.SiO2(i) <= 5 || isnan(selected_data.SiO2(i))) && ... % SiO2<=5
            (selected_data.TiO2(i) >= 30) && ... % 30<=TiO2
            (selected_data.FeO(i) >= 20) % 20<=FeO
        % Phase naming "Ilmenite"
        selected_data.Phase(i) = "Ilmenite";

        % "Magnetite"
        % 1) TiO2 <= 15 AND
        % 2) Al2O3 <= 5 AND
        % 3) 70 <= FeO
    elseif ismember('Magnetite',minerals) && ...
            (selected_data.TiO2(i) <= 15 || isnan(selected_data.TiO2(i))) && ... % TiO2<=15
            (selected_data.Al2O3(i) <= 5 || isnan(selected_data.Al2O3(i))) && ... % Al2O3<=5
            (selected_data.FeO(i) >= 70) % 70<=FeO
        % Phase naming "Magnetite"
        selected_data.Phase(i) = "Magnetite";

        % "Quartz"
        % 1) 90 <= SiO2
    elseif ismember('Quartz',minerals) && ...
        selected_data.SiO2(i) >= 90 % 90<=SiO2
        % Phase naming "Quartz"
        selected_data.Phase(i) = "Quartz";

        % "Spinel"
        % 1) SiO2 <= 5 AND
        % 2) 15 <= Al2O3 <= 60 AND
        % 3) CaO <= 1 AND
        % 4) 40 <= FeO+MgO+Cr2O3
    elseif ismember('Spinel',minerals) && ...
            (selected_data.SiO2(i) <= 5 || isnan(selected_data.SiO2(i))) && ... % SiO2<=5
            (selected_data.Al2O3(i) >= 15 && selected_data.Al2O3(i) <= 60) && ... % 15<=Al2O3<=60
            (selected_data.CaO(i) <= 1 || isnan(selected_data.CaO(i))) && ... % CaO<=1
            (sum([selected_data.FeO(i), selected_data.MgO(i), selected_data.Cr2O3(i)], 'omitnan') >= 40) % 40<=FeO+MgO+Cr2O3
        % Phase naming "Spinel"
        selected_data.Phase(i) = "Spinel";

        % "Plagioclase"
        % 1) 40 <= SiO2 <= 70 AND
        % 2) 15 <= Al2O3 <= 35 AND
        % 3) K2O <= 5 AND
        % 4) FeO+MnO+MgO <= 2 AMD
        % 5) 10 <= CaO+Na2O+K2O
    elseif ismember('Plagioclase',minerals) && ...
            (selected_data.SiO2(i) >= 40 && selected_data.SiO2(i) <= 70) && ... % 40<=SiO2<=70
            (selected_data.Al2O3(i) >= 15 && selected_data.Al2O3(i) <= 35) && ... % 15<=Al2O3<=35
            (selected_data.K2O(i) <= 5 || isnan(selected_data.K2O(i))) && ... % K2O<=5
            (sum([selected_data.FeO(i), selected_data.MnO(i), selected_data.MgO(i)], 'omitnan') <= 2) && ... % FeO+MnO++MgO<=2
            (sum([selected_data.CaO(i), selected_data.Na2O(i), selected_data.K2O(i)], 'omitnan') >= 10) % 10<=CaO+Na2O+K2O
        % Phase naming "Plagioclase"
        selected_data.Phase(i) = "Plagioclase";

        % "Mica"
        % 1) 30 <= SiO2 <= 55 AND
        % 2) 10 <= Al2O3 <= 30 AND
        % 3) 5 <= K2O
    elseif ismember('Mica',minerals) && ...
            (selected_data.SiO2(i) >= 30 && selected_data.SiO2(i) <= 55) && ... % 30<=SiO2<=55
            (selected_data.Al2O3(i) >= 10 && selected_data.Al2O3(i) <= 30) && ... % 10<=Al2O3<=30
            (selected_data.K2O(i) >= 5) % 5<=K2O
        % Phase naming "Mica"
        selected_data.Phase(i) = "Mica";

        % "Epidote"
        % 1) 30 <= SiO2 <= 45 AND
        % 2) 20 <= Al2O3 <= 30 AND
        % 3) 5 <= FeO <= 15 AND
        % 4) MnO <= 10 AND
        % 5) 12 <= CaO <= 25
    elseif ismember('Epidote',minerals) && ...
            (selected_data.SiO2(i) >= 30 && selected_data.SiO2(i) <= 45) && ... % 30<=SiO2<=45
            (selected_data.Al2O3(i) >= 20 && selected_data.Al2O3(i) <= 30) && ... % 20<=Al2O3<=30
            (selected_data.FeO(i) >= 5 && selected_data.FeO(i) <= 15) && ... % 5<=FeO<=15
            (selected_data.MnO(i) <= 10 || isnan(selected_data.MnO(i))) && ... % MnO<=10
            (selected_data.CaO(i) >= 12 && selected_data.CaO(i) <= 25) % 12<=CaO<=25
        % Phase naming "Epidote"
        selected_data.Phase(i) = "Epidote";

        % "Chlorite"
        % 1) 20 <= SiO2 <= 40 AND
        % 2) 15 <= Al2O3 <= 30 AND
        % 3) 10 <= FeO <= 25 AND
        % 4) 15 <= MgO <= 40
    elseif ismember('Chlorite',minerals) && ...
            (selected_data.SiO2(i) >= 20 && selected_data.SiO2(i) <= 40) && ... % 20<=SiO2<=40
            (selected_data.Al2O3(i) >= 10 && selected_data.Al2O3(i) <= 30) && ... % 15<=Al2O3<=30
            (selected_data.FeO(i) >= 10) && ... % 10<=FeO<=25
            (selected_data.MgO(i) >= 15 && selected_data.MgO(i) <= 40) % 15<=MgO<=40
        % Phase naming "Chlorite"
        selected_data.Phase(i) = "Chlorite";

        % "Olivine"
        % 1) 35 <= SiO2 <= 45 AND
        % 2) Al2O3 <= 5 AMD
        % 3) CaO <= 0.5 AND
        % 4) Na20 <= 0.5 AND
        % 5) K2O <= 0.5 AND
        % 6) 40 <= FeO+MgO
    elseif ismember('Olivine',minerals) && ...
            (selected_data.SiO2(i) >= 35 && selected_data.SiO2(i) <= 45) && ... % 35<=SiO2<=45
            (selected_data.Al2O3(i) <= 5) && ... % Al2O3<=5
            (selected_data.CaO(i) <= 0.5 || isnan(selected_data.CaO(i))) && ... % CaO<=0.5
            (selected_data.Na2O(i) <= 0.5 || isnan(selected_data.Na2O(i))) && ... % Na2O<=0.5
            (selected_data.K2O(i) <= 0.5 || isnan(selected_data.K2O(i))) && ... % K2O<=0.5
            (sum([selected_data.FeO(i), selected_data.MgO(i)], 'omitnan') >= 40) % 40<=FeO+MgO
        % Phase naming "Olivine"
        selected_data.Phase(i) = "Olivine";

        % "Cpx" - Clinopyroxene
        % 1) 45 <= SiO2 <= 60 AND
        % 2) Al2O3 <= 10 AND
        % 3) 10 <= MgO AND
        % 4) 5 <= CaO <= 25
        % 5) Na20 <= 1 AND
        % 6) K2O <= 0.5
    elseif ismember('Cpx',minerals) && ...
            (selected_data.SiO2(i) >= 45 && selected_data.SiO2(i) <= 60) && ... % 45<=SiO2<=60
            (selected_data.Al2O3(i) <= 10 || isnan(selected_data.Al2O3(i))) && ... % Al2O3<=10
            (selected_data.MgO(i) >= 10) && ... % 10<=MgO
            (selected_data.CaO(i) >= 5 && selected_data.CaO(i) <= 25) && ... % 5<=CaO<=25
            (selected_data.Na2O(i) <= 1 || isnan(selected_data.Na2O(i))) && ... % Na2O<=1
            (selected_data.K2O(i) <= 0.5 || isnan(selected_data.K2O(i))) % K2O<=0.5
        % Phase naming "Cpx"
        selected_data.Phase(i) = "Cpx";

        % "Opx" - Orthopyroxene
        % 1) 45 <= SiO2 <= 60 AND
        % 2) Al2O3 <= 10 AND
        % 3) 20 <= MgO AND
        % 4) CaO <= 10
        % 5) Na20 <= 0.5 AND
        % 6) K2O <= 0.5
        % 7) 30 <= FeO+MgO
    elseif ismember('Opx',minerals) && ...
            (selected_data.SiO2(i) >= 45 && selected_data.SiO2(i) <= 60) && ... % 45<=SiO2<=60
            (selected_data.Al2O3(i) <= 10 || isnan(selected_data.Al2O3(i))) && ... % Al2O3<=10
            (selected_data.MgO(i) >= 20) && ... % 20<=MgO
            (selected_data.CaO(i) <= 10 || isnan(selected_data.CaO(i))) && ... % CaO<=10
            (selected_data.Na2O(i) <= 0.5 || isnan(selected_data.Na2O(i))) && ... % Na2O<=0.5
            (selected_data.K2O(i) <= 0.5 || isnan(selected_data.K2O(i))) && ... % K2O<=0.5
            (sum([selected_data.FeO(i), selected_data.MgO(i)], 'omitnan') >= 30) % 30<=FeO+MgO
        % Phase naming "Opx"
        selected_data.Phase(i) = "Opx";

        % "Amphibole"
        % 1) 40 <= SiO2 <= 55 AND
        % 2) 1 <= Al2O3 <= 15 AND
        % 3) 5 <= CaO <= 15 AND
        % 4) 20 <= FeO+MgO
    elseif ismember('Amphibole',minerals) && ...
            (selected_data.SiO2(i) >= 40 && selected_data.SiO2(i) <= 55) && ... % 40<=SiO2<=55
            (selected_data.Al2O3(i) >= 1 && selected_data.Al2O3(i) <= 15) && ... % 1<=Al2O3<=15
            (selected_data.CaO(i) >= 5 && selected_data.CaO(i) <= 15) && ... % 5<=CaO<=15
            (sum([selected_data.FeO(i), selected_data.MgO(i)], 'omitnan') >= 20) % 20<=FeO+MgO
        % Phase naming "Amphibole"
        selected_data.Phase(i) = "Amphibole";

        % "Garnet"
        % 1) 35 <= SiO2 <= 45 AND
        % 2) 10 <= Al2O3 <= 25 AND
        % 3) 35 <= FeO+MnO+MgO+CaO
    elseif ismember('Garnet',minerals) && ...
            (selected_data.SiO2(i) >= 35 && selected_data.SiO2(i) <= 45) && ... % 35<=SiO2<=45
            (selected_data.Al2O3(i) >= 10 && selected_data.Al2O3(i) <= 25) && ... % 10<=Al2O3<=25
            (sum([selected_data.FeO(i), selected_data.MnO(i), selected_data.MgO(i), selected_data.CaO(i)], 'omitnan') >= 35) % 35<=FeO+MnO+MgO+CaO
        % Phase naming "Garnet"
        selected_data.Phase(i) = "Garnet";

        % If users want to add their own mineral phases,
        % it is recommended that they add them using "elseif" statements.

        % Data that does not meet the above conditions
        % is named "Undeterminable."
    else
        selected_data.Phase(i) = "Undeterminable";
    end
end

% Report identification mineral phases completion
disp('=== Mineral phases determination has been finished ===');

%% 10.Split the data (table) based on mineral phases
% Divide selected_data by mineral phases ("Phase").
% Use the local function (2) tableDivider for this process.
% If no data exists for a given mineral,
% an empty table is stored in the workspace.
% Export all the divided tables together as a single .xlsx file.
% NB: Save each mineral phases as a separate sheet.
% NB: No tables are exported if no data exists.
disp('------------------------------------------------------------------');
disp('=== Step 10: Dividing original table for each mineral phase ===');

% Devide the data table using local function (2) tableDivider
% Create a structure "table_data_struct",
% and store the data divided by mineral phases
table_data_struct = struct;
for i = 1:length(minerals)
    table_data_struct.(minerals{i}) = tableDivider(selected_data,selected_data.Phase,minerals{i});
end
% Store "Excluded" data in table_data_struct
table_data_struct.Excluded = tableDivider(selected_data,selected_data.Phase,"Excluded");
% Store "Undeterminable" data in table_data_struct
table_data_struct.Undeterminable = tableDivider(selected_data,selected_data.Phase,"Undeterminable");

% Export the split tables as a single .xlsx file
% Add the original (imported) filename to filename to export
filename_divided = [filename, '_Divided.xlsx'];
% Export all tables in table_data_struct
for i = 1:length(minerals)
    writetable(table_data_struct.(minerals{i}),filename_divided,'Sheet',minerals{i});
end
% Export the table of "Excluded" data in table_data_struct
writetable(table_data_struct.Excluded,filename_divided,'Sheet','Excluded');
% Export the table of "Undetermined" data in table_data_struct
writetable(table_data_struct.Undeterminable,filename_divided,'Sheet','Undetermined');

% Report exporting a file
disp(['Selected data written to: ', filename_divided]);

% Export dividing original table completion
disp('=== Dividing original table has been finished ===');

%% 11.Calculate elemental cation amounts
% Calculations are performed for each mineral phase (each separate table).
% Use the local function (3) calcCation.
% A separate table is assigned to each mineral phase.
% If no data exists for a particular mineral, 
% the table variable is not output.
disp('------------------------------------------------------------------');
disp('=== Step 11: Calculating element cations ===');

% Create a structure "table_cation_struct",
% and store the calculated cation data for each mineral phase here
table_cation_struct = struct;

% Calculate the amount of elemental cations 
% for all mineral phases stored in "minerals"
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        disp(['=== ' minerals{i} ' dataset ===']);
        table_cation_struct.(minerals{i}) = calcCation(table_data_struct.(minerals{i}),moduli_molarWeight,minerals{i});
    end
end

% Report calculating element cations completion
disp('=== Calculating element cations has been finished ===');

%% 12.Calculate the amount of moles of oxygen
% Calculations are performed for each mineral phase (each separate table).
% Use the local function (4) calcOxygen.
% A separate table is assigned to each mineral phase.
% If no data exists for a particular mineral, 
% the table variable is not output.
disp('------------------------------------------------------------------');
disp('=== Step 12: Calculating oxgen cations ===');

% Create a structure "table_cation_oxygen_struct",
% and store the calculated oxygen data for each mineral phase here
table_oxygen_struct = struct;

% Create a structure "oxegen_total_struct",
% and store the calculated oxygen total value for each mineral phase here
oxegen_total_struct = struct;

% Calculate the amount of oxygens [mol]
% for all mineral phases stored in "minerals"
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        disp(['=== ' minerals{i} ' dataset ===']);
        [table_oxygen_struct.(minerals{i}),oxegen_total_struct.(minerals{i})] = calcOxygen(table_data_struct.(minerals{i}),moduli_molarWeight);
    end
end

% Report calculating oxgen amounts completion
disp('=== Calculating oxgen amounts has been finished ===');

%% 13.Normalize cation data based on the assumed number of oxygen atoms
% Calculate the atomic per formula unit (a.p.f.u.).
% Calculations are performed for each mineral phase (each separate table).
% Use the local function (5) cationNormalize.
% If no data exists for a particular mineral, 
% the table variable is not output.
disp('------------------------------------------------------------------');
disp('=== Step 13: Normalizing based on the assumed number of oxygen atoms ===');

% Create a structure "table_cation_struct,
% and store the calculated oxygen data for each mineral phase here
table_cation_normalized_struct = struct;

% Calculate the a.p.f.u. values
% for all mineral phases stored in "minerals"
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        disp(['=== ' minerals{i} ' dataset ===']);
        table_cation_normalized_struct.(minerals{i}) = cationNormalize(table_cation_struct.(minerals{i}),oxegen_total_struct.(minerals{i}),moduli_molarWeight,moduli_Stoichiometry,minerals{i});
    end
end

% Report calculating a.p.f.u. completion
disp('=== Calculating a.p.f.u. has been finished ===');

%% 14.Output of cation calculation results
% The calculation results for cations, oxygen, and elements 
% are combined into a single table.
% A separate table is assigned to each mineral phase.
% If no data exists for a particular mineral, 
% the table variable is not output.
% The combined table is output as an .xlsx file.
disp('------------------------------------------------------------------');
disp('=== Step 14: Exporting the resuls of cation calculation ===');
disp('=== Combining tables and exporting to a single Excel sheet ===');


% Select the first column ("Comment") and combine them 
% for each mineral phase
% Create a structure "table_cation_oxygen_struct",
% and store the combined data for each mineral phase here
combined_table_struct = struct;
% Combine the tables (calculation results)
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        combined_table_struct.(minerals{i}) = [table_data_struct.(minerals{i}), table_cation_struct.(minerals{i})(:, 2:end), table_oxygen_struct.(minerals{i})(:, 2:end), table_cation_normalized_struct.(minerals{i})(:, 2:end)];
    end
end

% Export as an .xlsx file
% Add the original (imported) filename to filename to export
filename_cation_export = [filename, '_Results.xlsx'];

% Output one sheet per mineral phase
% Export all tables in combined_table_struct
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        writetable(combined_table_struct.(minerals{i}), filename_cation_export, 'Sheet', minerals{i});
    end
end

% Report exporting a file
disp(['Data has been exported to ' filename_cation_export]);

% Report exporting completion
disp('=== Exporting the resuls has been finished ===');

%% 15.Average values ​​for each sample code
% Average the calculation results 
% according to the sample code in the "Comment" column.
% Use the local function (14) columns2Mean.
% A separate table is assigned to each mineral phase.
% If no data exists for a particular mineral, 
% the table variable is not output.
% The combined table is output as an .xlsx file.
disp('------------------------------------------------------------------');
disp('=== Step 15: Averaging the values for every sample code ===');

% Extract data with the same sample code from the “Comment” column,
% and calculate the mean values
% Use the local function (14) columns2Mean
disp('=== Averaging all data for each mineral & sample code ===');

% Create a structure "table_data_mean_struct,
% and store the mean values for each sample code here
table_data_mean_struct = struct;

% Calculate the mean values
% for all mineral phases stored in "minerals"
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        table_data_mean_struct.(minerals{i}) = columns2Mean(table_data_struct.(minerals{i})); % measured data (wt%)
        table_cation_mean_struct.(minerals{i}) = columns2Mean(table_cation_struct.(minerals{i})); % calculated data (cations)
        table_cation_normalized_mean_struct.(minerals{i}) = columns2Mean(table_cation_normalized_struct.(minerals{i})); % normalized data (cations)
    end
end

% Summarize the average calculation results in one table
disp('=== Combining tables and exporting to a single Excel sheet ===');

% For each mineral phase, 
% combine all data except for the first column ("Comment")
% Create a structure "combined_table_mean_struct",
% and store the tables for each sample code here
combined_table_mean_struct = struct;

% Combine the averaged data for all mineral phases stored in "minerals"
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        combined_table_mean_struct.(minerals{i}) = [table_data_mean_struct.(minerals{i}), table_cation_mean_struct.(minerals{i})(:, 2:end), table_cation_normalized_mean_struct.(minerals{i})(:, 2:end)];
    end
end

% Export as an .xlsx file

% Add the original (imported) filename to filename to export
filename_cation_export = [filename, '_Mean.xlsx'];

% Output one sheet per mineral phase
% Export all tables in combined_table_struct
for i = 1:length(minerals)
    if ~isempty(table_data_struct.(minerals{i}))
        writetable(combined_table_mean_struct.(minerals{i}), filename_cation_export, 'Sheet', minerals{i});
    end
end

% Report exporting a file
disp(['Data has been exported to ' filename_cation_export]);

% Report exporting completion

disp('=== Averaging the result has been finished ===');

%% 16. Finish the script
% Report the end of scropt
disp('------------------------------------------------------------------');
disp('=== All processes/calculations have been finished ! ===');

% Acquire elapsed time
t_elapsed = toc;
fprintf('Processing time: %.4f 秒\n', t_elapsed)

end

%% ==================================
%% Local functions
%% ==================================
%% (1) apply_DL_mask
function masked_data = apply_DL_mask(data_matrix, DL_matrix)
% Syntax:
% function masked_data = apply_DL_mask(data_matrix, DL_matrix)
%
% Inputs:
% ------
% - data_matrix: @double
%     EPMA dataset for DL filtering
% - DL_matrix: @double
%     EPMA DL data for DL filtering
%
% Outputs:
% ------
% - masked_data
%     Filtered raw data based on DL
%
% Filter using Detection Limit (DL).
% Replace data below the DL value with NaN value.

% Check whether the number of input lines matches the number of DL lines
if size(data_matrix, 1) ~= size(DL_matrix, 1)
    % If they don't match, it will stop with an error
    error('Row number mismatch between selected_data and DL.');
end

% Convert values ​​below DL to NaN
masked_data = data_matrix;
mask = data_matrix < DL_matrix;
masked_data(mask) = NaN;

% Report removing the measured data completion
disp('=== All values less than DL were removed ===');

end

%% (2) tableDivider
function [table_divided] = tableDivider(table,column,target)
% Syntax:
% function [table_divided] = tableDivider(table,column,target)
%
% Inputs:
% ------
% - table: @double
%     EPMA dataset to be devided into sub-tables (based on mineral phases)
% - column:
%     criterion data within the table for table deviding
% - target:
%     standard data dor table deviding:
%     extract the values within the selected column that are the same as the target,
%     and output them as a new table.
%
% Outputs:
% ------
% - table_divided
%     devided table data due to column & target input
%
% Inputs a specific column, 
% groups together rows with the same data in that column, 
% and outputs them as a separate table.
% The "column" input is assumed to be the "Phase" column,
% and the "target" input is assumed to be the mineral phase
% (e.g., "Olivine").

% Extract the line number corresponding to "target"
target_id = ismember(column, target);

% create a table
table_divided = table(target_id,:);

% Report separating completion
disp(['=== ' target ' was separated into individual table ===']);

end

%% (3) calcCation
function [table_cation] = calcCation(table_data,moduli_molarWeight,mineralPhase)
% Syntax:
% function [table_cation] = calcCation(table_data,oxygen_total,moduli_molarWeight,mineralPhase)
%
% Inputs:
% ------
% - table_data: @table
%     EPMA dataset to be calculated into cations [mol]
% - oxegen_total @double
%     total amount of calculated oxygen [mol]
% - moduli_molarWeight: @cell
%     moduli of molarWeight for cation calculation: basically expected as Cation_moduli.xlsx
% - mineralPhase: @string
%     target mineral phase to be calculated into cations [mol]
%
% Outputs:
% ------
% - table_cation @table
%     table containing a result of calculation
%
% Calculate the amount (mol) of cations of each element
% from the measured data [wt%] of each oxide.
% For each mineral phase,
% assume the number of oxygen atoms in an ideal crystal lattice,
% and calculate the amount of cations in one crystal lattice.

% ------------------------------------------------------------------------
% 1. Setup
% Initialize a new table
table_cation = table();
% Add the first column ("Comment")
table_cation.Comment = table_data{:, 1};

% Check the location of "Total" column
id_total = find(strcmp("Total", table_data.Properties.VariableNames(:)));
tend = id_total-1;

% Check the size of "moduli_molarWeight"
csize = size(moduli_molarWeight);
cwidth = csize(2);

%Check the first column of "moduli_molarWeight" (Oxides)
% from the first column,
% and search for the same name as column_data
% Store the result (0 if not found)
match_indices = zeros(1,cwidth);

% Assign the variable name (element)
element_names = {};

% ------------------------------------------------------------------------
% 2. Calculate of cations & anions
% Calculate each column of "column_data" in turn
for col = 2:tend
    % Acquire the data in "table_data"
    % Column data (e.g., SiO2, TiO2)
    column_data = table_data{:, col};

    % Find an exact match one from the first line of "moduli_molarWeight"
    id = find(strcmp(table_data.Properties.VariableNames(col), moduli_molarWeight(1,:)));
    element = moduli_molarWeight{5,id};
    element_names{end+1} = element; %#ok<AGROW>

    % Record the matched location
    if ~isempty(id)
        match_indices(col) = id;
    else
        % If the corresponding column_data does not exist,
        % it will stop with an error
        error(['column_data does not exist: ', strjoin(column_data, ', ')]);
    end

    % Acquire cation_moduli data
    % Molar weights [g/mol]
    molar_data = moduli_molarWeight{2, id};

    % Number of cation atoms in the molecular formula of each oxide (Cation#)
    cation_number = moduli_molarWeight{3, id};

    % Check data types and calculate
    if isnumeric(column_data) && isnumeric(molar_data)
        % Create a variable to store the calculation result
        result = nan(size(column_data));
        % Index valid (non-NaN) values
        id_valid = ~isnan(column_data);
        % Calculate only the valid part
        % Formula: cation [mol] = measured data [wt%] / molar weight [g/mol] * number of cation
        result(id_valid) = column_data(id_valid) / molar_data * cation_number;
    else
        % If the data type is incorrect,
        % it will stop with an error
        error(['Invalid data type in column ' element]);
    end

    % Add results to the table
    table_cation.([element '_cation']) = result;
end

% ------------------------------------------------------------------------
% 3. Assign column names
% Obtain the variable name (column name) of "table_cation"
column_names = table_cation.Properties.VariableNames(:);

% For anions, use the suffix "_anion" instead of "_cation"
% Define a list of candidate anions (default: F and Cl)
anions = ["F_", "Cl_"];

% Match the variable names in column_names with "anions"
id_anion = contains(column_names, anions);

% Extract the part before the first "_" (element name)
column_names(id_anion) = strtok(column_names(id_anion), '_');

% Add the suffix "_anion" to the element name column name
column_names(id_anion) = strcat(column_names(id_anion), '_anion');
% Assign the column names
table_cation.Properties.VariableNames = column_names;

% ------------------------------------------------------------------------
% 4. Sum the calculated cations
% Calculate the total for each row and add it to the last column
disp("=== Summing rows for calculated cations & anions ===");

% Extract column names containing "_cation"
cation_vars = contains(table_cation.Properties.VariableNames, '_cation');

% Extract and sum only the data in the target column
% (row by row, ignoring NaNs)
row_sums_cation = sum(table_cation{:, cation_vars}, 2, 'omitnan');

% Add totals as a new column
table_cation.Total_cation = row_sums_cation;

% "Extract only columns with names containing "_anion"
anion_vars = contains(table_cation.Properties.VariableNames, '_anion');

% Extract and sum only the data in the target column
% (row by row, ignoring NaNs)
row_sums_anion = sum(table_cation{:, anion_vars}, 2, 'omitnan');

% Add totals as a new column
table_cation.Total_anion = row_sums_anion;

% ------------------------------------------------------------------------
% 5. Display the sums of cation
disp(['Cation sums of ' mineralPhase ':']);
disp(table_cation.Total_cation);

end

%% (4) calcOxygen
function [table_oxygen,oxygen_total] = calcOxygen(table_data,moduli_molarWeight)
% Syntax:
% function [table_cation_oxygen,oxygen_total] = calcOxygen(table_data,moduli_molarWeight)
%
% Inputs:
% ------
% - table_data: @table
%     EPMA dataset to be calculated into cations [mol]
% - moduli_molarWeight: @cell
%     moduli of molarWeight for cation calculation: basically expected as Cation_moduli.xlsx
%
% Outputs:
% ------
% - table_cation_oxygen @table
%     table containing a result of calculation
% - oxegen_total @double
%     total amount of calculated oxygen [mol]
%
% Calculate the amount of oxygen cations (mol)
% from the measured data [wt%] of each oxide.

% ------------------------------------------------------------------------
% 1. Setup
% Initialize a new table
table_oxygen = table();
% Add the first column ("Comment")
table_oxygen.Comment = table_data{:, 1};

% Check the location of "Total"
id_total = find(strcmp("Total", table_data.Properties.VariableNames(:)));
tend = id_total-1;

% Check the size of table_cation
csize = size(moduli_molarWeight);
cwidth = csize(2);

% Search the first column (Oxides) of table_cation,
% starting from the first, for a name with the same name as column_data
% Store the result (0 if not found)
match_indices = zeros(1,cwidth);

% Set the variable name (element)
element_names = {};

% ------------------------------------------------------------------------
% 2. Calculate oxygen molar
% cCalculate each column of "column_data" in turn
for col = 2:tend
    % Acquire column data of table_data (e.g., SiO2, TiO2)
    column_data = table_data{:, col};

    % Find an exact match from the first line of "moduli_molarWeight"
    id = find(strcmp(table_data.Properties.VariableNames(col), moduli_molarWeight(1,:)));
    element = moduli_molarWeight{5,id};

    % Element names are saved in order (the same order as table_data)
    element_names{end+1} = element; %#ok<AGROW>

    % Record the matched location
    if ~isempty(id)
        match_indices(col) = id;
    else
        % If the corresponding "column_data" does not exist,
        % it will stop with an error
        error(['column_data does not exist: ', strjoin(column_data, ', ')]);
    end

    % Obtain "cation_module" data
    % Molar weights [g/mol]
    molar_data = moduli_molarWeight{2, id};
    % The number of oxygen atoms in the molecular formula of each oxide (Oxigen#)
    oxygen_number = moduli_molarWeight{4, id};

    % Check data types and calculate (dealing with NaN)
    if isnumeric(column_data) && isnumeric(molar_data) && isnumeric(oxygen_number)
        % Create a variable to store the calculation result
        result = nan(size(column_data));
        % Index valid (non-NaN) values
        id_valid = ~isnan(column_data);
        % Calculate only the valid part
        % Formula: oxygen [mol] = measured data [wt%] / molar weight [g/mol]  * number of oxygen
        result(id_valid) = column_data(id_valid) / molar_data * oxygen_number;

    else
        % If the data type is incorrect,
        % it will stop with an error
        error(['Invalid data type in column ' element]);
    end

    % Add results to the table
    table_oxygen.([element '_oxygen']) = result;
end

% ------------------------------------------------------------------------
% 3. Sum the calculated oxygen
% Calculate the total for each row and add it to the last column
disp("=== Summing rows for calculated oxygen ===");

% Calculate row-by-row sums
row_sums = sum(table_oxygen{:, 2:end}, 2, 'omitnan');

% Add the total to the last column
table_oxygen.Total_oxygen = row_sums;

% ------------------------------------------------------------------------
% 4. Display the sums of cation
disp('Oxygen sums: ');
disp(table_oxygen.Total_oxygen);

% Total amount of oxygen cations (for output)
oxygen_total = table_oxygen.Total_oxygen;

end

%% (6) cationNormalize
function [table_cation_normalized] = cationNormalize(table_cation,oxygen_total,moduli_molarWeight,moduli_Stoichiometry,mineralPhase)
% Syntax:
% function [table_cation] = cationNormalize(table_data,oxygen_total,moduli_molarWeight,moduli_Stoichiometry,mineralPhase)
%
% Inputs:
% ------
% - table_data: @table
%     EPMA dataset to be calculated into cations [mol]
% - oxegen_total @double
%     total amount of calculated oxygen [mol]
% - moduli_molarWeight: @cell
%     moduli of molarWeight for cation calculation: basically expected as Cation_moduli.xlsx
% - moduli_Stoichiometry: @table
%     moduli of stoichiometry for cation calculation: basically expected as Cation_moduli.xlsx
% - mineralPhase: @string
%     target mineral phase to be calculated into cations [mol]
%
% Outputs:
% ------
% - table_cation @table
%     table containing a result of calculation
%
% Calculate the amount (mol) of cations of each element
% from the measured data [wt%] of each oxide.
% For each mineral phase,
% assume the number of oxygen atoms in an ideal crystal lattice,
% and calculate the amount of cations in one crystal lattice.

% ------------------------------------------------------------------------
% 1. Setup
% Extract the number of oxygen atoms in one crystal lattice (Oxygen#)
% for each mineral phase as a variable.
oxygen_number = moduli_Stoichiometry.(mineralPhase)(1);

% Initialize a new table
table_cation_normalized = table();
% Add the first column ("Comment")
table_cation_normalized.Comment = table_cation{:, 1};

% Check the location of "Total"
id_total = find(strcmp("Total_cation", table_cation.Properties.VariableNames(:)));
tend = id_total-1;

% Check the size of moduli_molarWeight
csize = size(moduli_molarWeight);
cwidth = csize(2);

% Change the name of the fifth line of moduli_molarWeight
moduli_molarWeight_names = moduli_molarWeight(5,:);
moduli_molarWeight_names = strcat(moduli_molarWeight_names, '_cation');

% Check the first column (Oxides) of moduli_molarWeight
% from the first column and search for the same name as column_data
% Store the result (0 if not found)
match_indices = zeros(1,cwidth);

% Assign the variable name (element)
cation_names = {};

% ------------------------------------------------------------------------
% 2. Normalize cations
% The ratio of elemental cations to the original oxide varies
% (depending on the element)
% ==> This can be handled by the module
for col = 2:tend
    % Acquire the data for table_cation
    column_data = table_cation{:, col};  % 列データ（e.g., Si_cation, Ti_cation)

    % Find an exact match from the first line of moduli_molarWeight
    id = find(strcmp(table_cation.Properties.VariableNames(col), moduli_molarWeight_names(:)));
    if id ~= 0
        cation_name = moduli_molarWeight_names{id};


        cation_names{end+1} = cation_name; %#ok<AGROW>

        % Record the matched location
        if ~isempty(id)
            match_indices(col) = id;
        else
            % If the corresponding "column_data" does not exist,
            % it will stop with an error
            error(['column_data does not exist: ', strjoin(column_data, ', ')]);
        end

        % Check data types and calculate
        if isnumeric(column_data) && isnumeric(oxygen_number)
            % Create a variable to store the calculation result
            result = nan(size(column_data));
            % Index valid (non-NaN) values
            id_valid = ~isnan(column_data) & ~isnan(oxygen_total);
            % Calculate only the valid part
            % Formula: cation_normalized [mol] = calculated cations [mol] * number of oxygen atoms [mol]  / total oxygen molar [mol]
            result(id_valid) = column_data(id_valid) ./ oxygen_total(id_valid) * oxygen_number;
        else
            % If the data type is incorrect,
            % it will stop with an error
            error(['Invalid data type in column ' num2str(col)]);
        end

        % Add results to the table
        table_cation_normalized.([cation_name '_normalized']) = result;
    end
end

% ------------------------------------------------------------------------
% 3. Calculate Fe3+ cations
% Use local function (6) calcFe3

% For minerals required Fe3+ calculation
if ismember(string(mineralPhase), ["Spinel", "Magnetite", "Opx", "Cpx", "Epidote", "Garnet"])
    [Fe3, Fe2_calcFe3] = calcFe3(table_cation_normalized, moduli_molarWeight, mineralPhase);
    table_cation_normalized.Fe_cation_normalized = Fe2_calcFe3;
else
    % For minerals NOT required Fe3+ calculation
    % Check the size of "table_cation"
    h = size(table_cation_normalized);
    % Create a dummy column using the height of the table_cation (h(1)).
    Fe3 = NaN(h(1),1);
end

% Add Fe3+ & Fe2 values to the table
table_cation_normalized.Fe3_cation_normalized = Fe3;

% ------------------------------------------------------------------------
% 4. Sum the calculated cations
% Calculate the total for each row and add it to the last column
disp("=== Summing rows for normalized cations ===");

% Extract column names containing "_cation"
cation_vars = contains(table_cation_normalized.Properties.VariableNames, '_normalized');

% Extract and sum only the data in the target column
% (row by row, ignoring NaNs)
row_sums_cation = sum(table_cation_normalized{:, cation_vars}, 2, 'omitnan');

% Add totals as a new column
table_cation_normalized.Total_cation_normalized = row_sums_cation;

% ------------------------------------------------------------------------
% 5. Displaying the sums of cation
disp(['Normalized cation sums of ' mineralPhase ':']);
disp(table_cation_normalized.Total_cation_normalized);

% ------------------------------------------------------------------------
% 6. Calculate petrological factors
% Calculate of Mg#, Cr#, F value, An content, etc.
% Olivine
if mineralPhase == "Olivine"
    % use local function (7)
    table_cation_normalized = calcOl(table_cation_normalized);

    % Opx
elseif mineralPhase == "Opx"
    % use local function (8)
    table_cation_normalized = calcOpx(table_cation_normalized);

    % Cpx
elseif mineralPhase == "Cpx"
    % use local function (9)
    table_cation_normalized = calcCpx(table_cation_normalized);

    % Spinel
elseif mineralPhase == "Spinel"
    % use local function (10)
    table_cation_normalized = calcSp(table_cation_normalized);

    % Plagioclase
elseif mineralPhase == "Plagioclase"
    % use local function (11)
    table_cation_normalized = calcPl(table_cation_normalized);

    % Amphibole
elseif mineralPhase == "Amphibole"
    % use local function (12)
    table_cation_normalized = calcAmph(table_cation_normalized);

    % Garnet
elseif mineralPhase == "Garnet"
    % use local function (13)
    table_cation_normalized = calcGrt(table_cation_normalized);
end

end

%% (6) calcFe3
function [Fe3, Fe2] = calcFe3(table_cation, moduli_molarWeight, mineralPhase)
% Syntax:
% function [Fe3, Fe2] = calcFe3(table_cation, moduli_molarWeight, mineralPhase)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
% - moduli_molarWeight: @cell
%     moduli of molarWeight for cation calculation: basically expected as Cation_moduli.xlsx
% - mineralPhase: @string
%     target mineral phase to be calculated into cations [mol]
%
% Outputs:
% ------
% - Fe3 @double
%     molar amount of Fe3+
% - Fe2 @double
%     molar amount of Fe2+
%
% Calculate Fe3+ from FeO (Fe2+).
% Calculate separately for each mineral phase

% Find the location of FeO and Fe2O3 from "moduli_molarWeight"
id_feo = strcmp("FeO", moduli_molarWeight(1,:));
id_fe2o3 = strcmp("Fe2O3", moduli_molarWeight(1,:));

% Extract the molar weights of FeO and Fe2O3 from "moduli_molarWeight"
mw_feo = cell2mat(moduli_molarWeight(2,id_feo));
mw_fe2o3 = cell2mat(moduli_molarWeight(2,id_fe2o3));

% ------------------------------------------------------------------------
% "MgAl2O4" - Fe3O4 series
% Droop (1987)
if ismember(string(mineralPhase), ["Spinel", "Magnetite"])
    % Calculate Fe3+
    table_cation = fillmissing(table_cation(:,2:end), 'constant', 0);
    A = 2 * sum(fillmissing(table_cation{:, {'Fe_cation_normalized', 'Mn_cation_normalized','Mg_cation_normalized','Ca_cation_normalized','Na_cation_normalized','K_cation_normalized', 'Ni_cation_normalized'}}, 'constant', 0), 2);
    B = sum(fillmissing(table_cation{:, {'Si_cation_normalized','Ti_cation_normalized','Al_cation_normalized','Cr_cation_normalized','V_cation_normalized'}}, 'constant', 0), 2);
    Fe3 =  (A - B)/3;

    % If the calculated Fe3+ shows a negative value,
    % the value is considered to be 0
    id_fe3minus = ismember(Fe3, Fe3(Fe3<0));
    Fe3(id_fe3minus) = 0;

    % Calculate Fe2+
    Fe2 = fillmissing(table_cation.Fe_cation_normalized, 'constant', 0) - Fe3;

    % If the calculated Fe3+ shows a negative value,
    % the value is considered to be 0
    id_fe2minus = ismember(Fe2, Fe2(Fe2<0));
    Fe2(id_fe2minus) = 0;

    % Opptional: If assuming Fe3+ = 0
    % % Calculate Fe3+
    % Fe3 = zeros(height(table_cation(:,1)),1);
    % 
    % % Calculate Fe2+
    % Fe2 = table_cation.Fe_cation;

    % ------------------------------------------------------------------------
    % "Pyroxenes"
    % Papike et al. (1974)
elseif ismember(string(mineralPhase), ["Opx", "Cpx"])
    % Calculate Fe3+
    SiteCharge_Deficiency = sum(fillmissing(table_cation{:, {'Al_cation_normalized','Na_cation_normalized'}}, 'constant', 0), 2);
    SiteCharge_Excess = sum(fillmissing(table_cation{:, {'Al_cation_normalized','Cr_cation_normalized','Ti_cation_normalized','Ti_cation_normalized'}}, 'constant', 0), 2);
    Fe3 = SiteCharge_Deficiency - SiteCharge_Excess;

    % If the calculated Fe3+ shows a negative value,
    % the value is considered to be 0
    id_fe3minus = ismember(Fe3, Fe3(Fe3<0));
    Fe3(id_fe3minus) = 0;

    % Calculate Fe2+
    Fe2 = fillmissing(table_cation.Fe_cation_normalized, 'constant', 0) - Fe3;

    % If the calculated Fe2+ shows a negative value,
    % the value is considered to be 0
    id_fe2minus = ismember(Fe2, Fe2(Fe2<0));
    Fe2(id_fe2minus) = 0;

    % Opptional: If assuming Fe3+ = 0
    % % % Calculate Fe3+
    % % Fe3 = zeros(height(table_cation(:,1)),1);
    % % 
    % % % Calculate Fe2+
    % % Fe2 = table_cation.Fe_cation;

    % ------------------------------------------------------------------------
    % "Garnet"
    % Enami (2012)
elseif ismember(string(mineralPhase), "Garnet")
    % Calculate Fe3+
    A = 2 * sum(fillmissing(table_cation{:, {'Si_cation_normalized','Ti_cation_normalized', 'Fe_cation_normalized', 'Mn_cation_normalized','Mg_cation_normalized','Ca_cation_normalized', 'Ni_cation_normalized'}}, 'constant', 0), 2);
    B = 3 * sum(fillmissing(table_cation{:, {'Al_cation_normalized','Cr_cation_normalized','Na_cation_normalized','K_cation_normalized','V_cation_normalized'}}, 'constant', 0), 2);
    Fe3 = (A - B) / 5;

    % If the calculated Fe3+ shows a negative value,
    % the value is considered to be 0
    id_fe3minus = ismember(Fe3, Fe3(Fe3<0));
    Fe3(id_fe3minus) = 0;

    % Calculate Fe2+
    Fe2 = fillmissing(table_cation.Fe_cation_normalized, 'constant', 0) - Fe3;

    % If the calculated Fe2+ shows a negative value,
    % the value is considered to be 0
    id_fe2minus = ismember(Fe2, Fe2(Fe2<0));
    Fe2(id_fe2minus) = 0;

    % Opptional: If assuming Fe3+ = 0
    % % Calculate Fe3+
    % Fe3 = zeros(height(table_cation(:,1)),1);
    % 
    % % Calculate Fe2+
    % Fe2 = table_cation.Fe_cation_normalized;

    % ------------------------------------------------------------------------
    % "Epidote"
    % Masumoto et al. (2014): Total iron was treated as Fe2O3
elseif ismember(string(mineralPhase), "Epidote")
    % Calculate Fe3+
    Fe3 = table_cation.Fe_cation_normalized;

    % Calculate Fe2+
    Fe2 = zeros(height(table_cation(:,1)),1);

    % Warning message for assuming Fe2+ >>> Fe3+
    warning([mineralPhase ' is  intended to be Fe3+ >>> Fe2+'])

    % Opptional: If assuming Fe3+ = 0
    % % Calculate Fe3+
    % Fe3 = zeros(height(table_cation(:,1)),1);
    % 
    % % Calculate Fe2+
    % Fe2 = table_cation.Fe_cation_normalized;

    % ------------------------------------------------------------------------
else
    % Calculate Fe3+
    Fe3 = zeros(height(table_cation(:,1)),1);

    % Calculate Fe2+
    Fe2 = table_cation.Fe_cation;

    % Warning message for Fe3+ calculation not expected
    warning([mineralPhase ' is not intended to calculate Fe3+'])

end

% Report Fe3+ calculation completion
disp(['=== Fe3+ calculation for ' mineralPhase ' has been finished ===']);

end

%% (7) calcOl
function table_cation = calcOl(table_cation)
% Syntax:
% function table_cation = calcOl(table_cation)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Outputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Calculate various parameters of Olivine
% Store the calculation results in "table_cation" and output them

% Mg# for olivine
% Mg# = Mg/(Mg+Fe2+)
Mgnum_ol = table_cation.Mg_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Mg#') = Mgnum_ol;

% Fe# for olivine
% Fe# = Fe/(Mg+Fe2+)
Fenum_ol = table_cation.Fe_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Fe#') = Fenum_ol;

end

%% (8) calcOpx
function table_cation = calcOpx(table_cation)
% Syntax:
% function table_cation = calcOpx(table_cation)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Outputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Calculate various parameters of Opx
% Store the calculation results in "table_cation" and output them

% Mg# for opx
% Mg# = Mg/(Mg+Fe2+)
Mgnum_opx = table_cation.Mg_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Mg#') = Mgnum_opx;

% Fe# for opx
% Fe# = Fe/(Mg+Fe2+)
Fenum_opx = table_cation.Fe_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Fe#') = Fenum_opx;

end

%% (9) calcCpx
function table_cation = calcCpx(table_cation)
% Syntax:
% function table_cation = calcCpx(table_cation)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Outputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Calculate various parameters of Cpx
% Store the calculation results in "table_cation" and output them

% Mg# for cpx
% Mg# = Mg/(Mg+Fe2+)
Mgnum_cpx = table_cation.Mg_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Mg#') = Mgnum_cpx;

% Fe# for cpx
% Fe# = Fe/(Mg+Fe2+)
Fenum_cpx = table_cation.Fe_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Fe#') = Fenum_cpx;

end

%% (10) calcSp
function table_cation = calcSp(table_cation)
% Syntax:
% function table_cation = calcSp(table_cation)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Outputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Calculate various parameters of Spinel
% Store the calculation results in "table_cation" and output them

% Mg# for spinel
% Mg# = Mg/(Mg+Fe2+)
Mgnum_sp = table_cation.Mg_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Mg#') = Mgnum_sp;

% Fe# for spinel
% Fe# = Fe2+/(Mg+Fe2+)
Fenum_sp = table_cation.Fe_cation_normalized./(table_cation.Mg_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.('Fe#') = Fenum_sp;

% Cr# for spinel
% Cr# = Cr/(Cr+Al)
Crnum_sp = table_cation.Cr_cation_normalized./(table_cation.Cr_cation_normalized + table_cation.Al_cation_normalized);
table_cation.('Cr#') = Crnum_sp;

% YCr for spinel
% YCr = Cr/(Cr+Al+Fe3+)
YCr = table_cation.Cr_cation_normalized./(table_cation.Cr_cation_normalized + table_cation.Al_cation_normalized + table_cation.Fe3_cation_normalized);
table_cation.YCr = YCr;

% YFe for spinel
% YFe = Fe3+/(Cr+Al+Fe3+)
YFe = table_cation.Fe3_cation_normalized./(table_cation.Cr_cation_normalized + table_cation.Al_cation_normalized + table_cation.Fe3_cation_normalized);
table_cation.YFe = YFe;

% RFeCr
% RFeCr =-- Fe3+/(Fe3++Cr)
RFeCr = table_cation.Fe3_cation_normalized./(table_cation.Fe3_cation_normalized + table_cation.Cr_cation_normalized);
table_cation.RFeCr = RFeCr;

% RFe
% RFe = Fe3+/(Fe3++Fe2+)
RFe = table_cation.Fe3_cation_normalized./(table_cation.Fe3_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.RFe = RFe;

% RFe3+Fe2+
% RFe3+2+ = Fe3+/Fe2+
RFe3Fe2 = table_cation.Fe3_cation_normalized./table_cation.Fe_cation_normalized;
table_cation.RFe3Fe2 = RFe3Fe2;

% Fvalue (degree of melting) by Hellebrand et al. (2001)
% F = 10*ln(Cr#) + 24
Fvalue = 10*log(Crnum_sp) + 24;
table_cation.Fvalue = Fvalue;

end

%% (11) calcPl
function table_cation = calcPl(table_cation)
% Syntax:
% function table_cation = calcPl(table_cation)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Outputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Calculate various parameters of Plagioclase
% Store the calculation results in "table_cation" and output them

% Ca, Na, Kの合計値
total_AnAbOr = table_cation.Ca_cation_normalized + table_cation.Na_cation_normalized + table_cation.K_cation_normalized;

% An% for plagioclase
% An% = Ca/(Ca+Na+K)
An_pl = table_cation.Ca_cation_normalized./total_AnAbOr;
table_cation.An = An_pl;

% Ab% for plagioclase
% Ab% = Na/(Ca+Na+K)
Ab_pl = table_cation.Na_cation_normalized./total_AnAbOr;
table_cation.Ab = Ab_pl;

% Or% for plagioclase
% Or% = K/(Ca+Na+K)
Or_pl = table_cation.K_cation_normalized./total_AnAbOr;
table_cation.Or = Or_pl;

end

%% (12) calcAmph
function table_cation = calcAmph(table_cation)
% Syntax:
% function table_cation = calcAmph(table_cation)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Outputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Calculate various parameters of Amphibole
% Store the calculation results in "table_cation" and output them
%
% ------------------------------------------------------------------------
% Nomenclature of Amphiboles:
% LEAKE, B.E. et al. (1997), The Canadian Mineralogist Vol.35, 219-246.
%
% standard amphibole formula: A B2(VI) C5(IV) T8 O22 (OH)2
%
% ------------------------------------------------------------------------
% The components of the formula conventionally described as
% A, B, C, T and “OH” correspond to the following crystallographic sites:
%
% A: one site per formula unit
% B: two M4 sites per formula unit
% C: a composite of five sites made up of two M1, two M2 & one M3 sites
%    per formula unit
% T: eight sites, in two sets of four, which need not be distinguished
%    in this document;
%  “OH”: two sites per formula unit
%
% ------------------------------------------------------------------------
% The ions considered "normally" to occupy these sites
% are in the following categories:
%
% [](empty site) & K: @ A site only
% Na: @ A or B site
% Ca: @ B site only
% L-type ions: @ C or B site
%   *Mg, Fe2+, Mn2+, Li, & rarer ions of similar size, such as Zn, Ni, Co
% M-type ions: @ C site only
%   *Al, Fe3+, & more rarely, Mn3+,Cr3+
%   **M-type ions normally occupy M2 sites;
%     so are normally limited to two of the five C sites
% High-valency ions: @ C or T site
%   *Ti4+
% Zr4+: @ C site only
% Si: @ T site only
% Anions: @ "OH" site
%   *OH, F, CL, O
%
% ------------------------------------------------------------------------
% Calculation flow for standard amphibole formula
% recommended by Leake et al. (1997)
% *NB: In this calculation, uncertainty results from
%      lack of analyses for H2O, Fe3+ and Fe2+:
%
% 1) If H2O and halogen contents are well established:
%    the formula should be calculated to 24(O,OH,F,Cl)
%    *NB: this script assumed 23(O) following 2) rule
% 2) If the H2O plus halogen content is uncertain:
%    the formula should be calculated to the basis of 23(O)
%    with 2(OH,F,Cl) assumed
%    *unless this leads to an impossibility of
%     satisfying any of the  following criteria,
%     in which case an appropriate change in the assumed number of
%     (OH + F + Cl) should be made
%    **In this script, 23(O) & 0(OH + F + Cl) are assumed
% 3) Sum T site to 8.00 using Si, then Al, then Ti
%    *Fe3+ is not allocated to T site
%    **The normal maximum substitution for Si is 2, but it can be exceeded
% 4) Sum C site to 5.00 using excess Al and Ti from (3)
%    and then successively Zr, Cr3+,Fe3+,Mn3+, Mg, Fe2+, Mn2+,
%    any other L2+-type ions, and then Li
% 5) Sum B site to 2.00 using excess Mg, Fe2+,Mn2+ and Li from (4)
%    then Ca, then Na
% 6) Excess Na from (5) is assigned to A site, then all K
%    *Total A site should be between 0 and 1.00
%
% ------------------------------------------------------------------------
% Four amphibole  classification depending on the occupancy  of the B sites
%
% 1) magnesium - iron - manganese - lithium group
%    (Ca + Na)_B < 1.00 && the sum of L-type ions (Mg,Fe,Mn,Li)_B => 1.00
% 2) calcic group
%    (Ca + Na)_B => 1.00 && Na_B < 0.50
%    *Usually, but not in every case, Ca_B > 1.50
% 3) sodic–calcic group
%    (Ca + Na)B > 1.00 && 0.50 =< Na_B < 1.50
% 4) sodic group (previously referred to alkali amphiboles)
%    Na_B => 1.50
%
% ------------------------------------------------------------------------
%
% Check the size of table_cation
stable = size(table_cation);
w = stable(2);

% Calculate occupancy rate, etc.
% Group A
% Site occupancy across amphibole
f1 = 16 ./ table_cation.Total_cation_normalized;
% Si occupancy at T site
f2 = 8 ./ table_cation.Si_cation_normalized;
% Site occupancy at T, B, and C sites (excluding A site)
f3 = 15 ./ (table_cation.Total_cation_normalized - table_cation.Na_cation_normalized - table_cation.K_cation_normalized + table_cation.Mn_cation_normalized);
% Ca occupancy at the B site
f4 = 2 ./ table_cation.Ca_cation_normalized;
% Lower limit check
f5 = ones(length(f1), 1);

% Group B
% Site occupancy of Si and Al at the T site
f6 = 8 ./ (table_cation.Si_cation_normalized + table_cation.Al_cation_normalized);
% Site occupancy rates at T, B, and C sites (excluding K at A site)
f7 = 15 ./ (table_cation.Total_cation_normalized - table_cation.K_cation_normalized);
% Site occupancy rate at T and C sites (excluding A and B sites)
f8 = 13 ./ (table_cation.Total_cation_normalized - table_cation.Ca_cation_normalized - table_cation.Na_cation_normalized - table_cation.K_cation_normalized);
% 10Fe³⁺ Scaling factor for normalization N
f9 = 36 ./ (46 - table_cation.Si_cation_normalized - table_cation.Ti_cation_normalized - table_cation.Al_cation_normalized - table_cation.Cr_cation_normalized);
% Contribution of Fe
f10 = 46 ./ (table_cation.Fe_cation_normalized + 46);

% Check integrity
fmins = [];
fmaxs = [];
for i = 1:length(f1)
    fmin = min([f1(i) f2(i) f3(i) f4(i) f5(i)]);
    fmins = [fmins; fmin]; %#ok<AGROW>
    fmax = max([f6(i) f7(i) f8(i) f9(i) f10(i)]);
    fmaxs = [fmaxs;fmax ]; %#ok<AGROW>
    % There is no problem if fmax < 1
    % NB: An occupancy greater than 1 is physically inappropriate
    % (excessive cation allocation)
end

% Average of fmin and fmax
meanAB = (fmins + fmaxs) / 2;

% Uniform rescaling with a correction factors
for i = 1:w
    if contains(table_cation.Properties.VariableNames{i}, '_cation_normalized')
        % Multiply the cation data by the correction coefficient "meanAB"
        cation4Amph = table_cation{:,i} .* meanAB;
        % Save results in the table as element name + _cation4Amph
        cation4Amph_name = char(table_cation.Properties.VariableNames{i});
        cation4Amph_name = strtok(cation4Amph_name, '_');
        table_cation.([cation4Amph_name '_cation4Amph']) = cation4Amph;
    end
end

% Estimate Fe3+
Fe3_cation4Amph = 46 * (1 - meanAB);
table_cation.Fe3_cation4Amph = Fe3_cation4Amph;
% Calculate Fe2+
Fe2_cation4Amph = table_cation.Fe_cation4Amph - Fe3_cation4Amph;
table_cation.Fe2_cation4Amph = Fe2_cation4Amph;

% cation minor
cm = table_cation.Total_cation4Amph - table_cation.Ca_cation4Amph- table_cation.Na_cation4Amph - table_cation.K_cation4Amph - 13;

% Calculate the moles per site for each element
% T site (Tetrahedral cations)
% Si at T site
Si_T = table_cation.Si_cation4Amph;
table_cation.Si_T = Si_T;
% Al(IV) at T site
Al_T = 8 - table_cation.Si_cation4Amph;
table_cation.Al_T = Al_T;
% Total amount of T sites
% NB: Calculate the total to 8
Total_T = Si_T + Al_T;
table_cation.Total_T = Total_T;

% C site
% Al(IV) at C site
Al_C = table_cation.Al_cation4Amph - Al_T;
table_cation.Al_C = Al_C;
% Ti at C site
Ti_C = table_cation.Ti_cation4Amph;
table_cation.Ti_C = Ti_C;
% Fe3+ at C site
Fe3_C = table_cation.Fe3_cation4Amph;
table_cation.Fe3_C = Fe3_C;
% Cr at C site
Cr_C = table_cation.Cr_cation4Amph;
table_cation.Cr_C = Cr_C;
% Ni at C site
Ni_C = table_cation.Ni_cation4Amph;
table_cation.Ni_C = Ni_C;
% V at C site
V_C = table_cation.V_cation4Amph;
table_cation.V_C = V_C;
% Mg at C site
Mg_C = table_cation.Mg_cation4Amph;
table_cation.Mg_C = Mg_C;
% Fe2+ at C site
Fe2_C = 5 - (Al_C + Ti_C + Fe3_C + Cr_C + Ni_C + V_C + Mg_C);
table_cation.Fe2_C = Fe2_C;
% Mn at C site
Mn_C = 5 - (Al_C + Ti_C + Fe3_C + Cr_C + Ni_C + V_C + Mg_C + Fe2_C);
table_cation.Mn_C = Mn_C;
% Total amount of C sites
Total_C = Al_C + Ti_C + Fe3_C + Cr_C + Ni_C + V_C + Mg_C + Fe2_C + Mn_C;
table_cation.Total_C = Total_C;

% B Site
% Fe2+ at B site (if the value is negative, it will be set to 0.)
Fe2_B = table_cation.Fe2_cation4Amph - Fe2_C;
if Fe2_B < 0
    Fe2_B = 0;
end
table_cation.Fe2_B = Fe2_B;
% Mn at B site
Mn_B = table_cation.Mn_cation4Amph - Mn_C;
table_cation.Mn_B = Mn_B;
% Ca at B site
Ca_B = table_cation.Ca_cation4Amph;
table_cation.Ca_B = Ca_B;
% Na at B site
Na_B = 2 - (Fe2_B + Mn_B + Ca_B);
table_cation.Na_B = Na_B;
% Total amount of B sites
Total_B = Fe2_B + Mn_B + Ca_B + Na_B;
table_cation.Total_B = Total_B;
% Ca & Na at B site
CaNa_B = Ca_B + Na_B;
table_cation.CaNa_B = CaNa_B;

% A Site
% Na at A site
% NB: Assuming that the B site (M4 site) is completely filled with two cations
Na_A = table_cation.Na_cation4Amph - Na_B;
table_cation.Na_A = Na_A;
% K at A site
K_A = table_cation.K_cation4Amph;
table_cation.K_A = K_A;
% Na & K at A site (Total amount of A sites)
NaK_A = Na_A + K_A;
table_cation.NaK_A = NaK_A;

% Total amount of cations after site allocation
Total_sites = Total_T + Total_C + Total_B + NaK_A;
table_cation.Total_sites = Total_sites;


% Detailed site allocation
% T1 & T2 site Allocation
% Al(IV) site allocation
% Estimate Si in the T1 site
% NB: Si preferentially coordinates to the T1 site
Si_T1 = table_cation.Si_cation4Amph - 4;
table_cation.Si_T1 = Si_T1;
% EstimateAl(IV) at T1 site
Al_T1 = 4 - Si_T1;
table_cation.Al_T1 = Al_T1;
% Estimate remaining Al(IV)
% NB: presumably in the M2 site in the C site
Al_M2 = table_cation.Al_cation4Amph - Al_T1;
table_cation.Al_M2 = Al_M2;

% cation minor at M4 site
cm4 = cm;

% M4 site Allocation
% Ca at M4 site
Ca_M4 = table_cation.Ca_cation4Amph;
table_cation.Ca_M4 = Ca_M4;
% Na at M4 site
Na_M4 = table_cation.Na_cation4Amph + table_cation.K_cation4Amph +15 - table_cation.Total_cation4Amph;
table_cation.Na_M4 = Na_M4;

% A site Allocation
% Na at A site
% NB: Assuming that all 15 T+C+B sites are filled
Na_A2 = table_cation.Total_cation4Amph - table_cation.K_cation4Amph - 15;
table_cation.Na_A2 = Na_A2;
% K at A site
K_A2 = table_cation.K_cation4Amph;
table_cation.K_A2 = K_A2;
% Vacant at A site
vA = 16 - table_cation.Total_cation4Amph;
table_cation.vA = vA;

% Calculate Mg#
Mgnum = table_cation.Mg_cation4Amph ./ (table_cation.Mg_cation4Amph + table_cation.Fe2_cation4Amph);
table_cation.Mgnum = Mgnum;

% M1 & M3 site Allocation
% Calculate the amount of Mg at C site (M1+M3 site) based on Mg#
% Assume Mg + Fe²⁺ = 3 (M1 + M2)
Mg_M1M3 = 3 * Mgnum;
table_cation.Mg_M1M3 = Mg_M1M3;
% Calculate the amount of Fe2+ at C site (M1+M3 site)
% Assuming that Fe2+ is allocated to the C site in a complementary manner with Mg
Fe2_M1M3 = 3 * Mg_M1M3;
table_cation.Fe2_M1M3 = Fe2_M1M3;

% M2 site Allocation
% Fe3+ at M2 site
Fe3_M2 = Fe3_cation4Amph;
table_cation.Fe3_M2 = Fe3_M2;
% Ti at M2 site
Ti_M2 = Ti_C;
table_cation.Ti_M2 = Ti_M2;
% Calculate the amount of Mg at C site (M2 site) based on Mg#
% 10 is the total number of cations in the T site + C site + B site
% (or the sum of the M1 to M3 sites)
Mg_M2 = Mgnum .* (10 - table_cation.Si_cation4Amph - table_cation.Ti_cation4Amph - table_cation.Al_cation4Amph - table_cation.Fe3_cation4Amph);
table_cation.Mg_M2 = Mg_M2;
% Calculate the amount of Fe2+ at C site (M2 site)
% Assuming that Fe2+ is allocated to the C site in a complementary manner with Mg
Fe2_M2 = (1 - Mgnum) .* (10 - table_cation.Si_cation4Amph - table_cation.Ti_cation4Amph - table_cation.Al_cation4Amph - table_cation.Fe3_cation4Amph);
table_cation.Fe_M2 = Fe2_M2;

% Total cations at each site
% M4 site
Total_M4 = cm4 + Ca_M4 + Na_M4;
table_cation.Total_M4 = Total_M4;
% A site
Total_A = vA + K_A2 + Na_A2;
table_cation.Total_A = Total_A;
% T1 site
Total_T1 = Si_T1 + Al_T1;
table_cation.Total_T1 = Total_T1;
% M2 site
Total_M2 = Al_M2 + Fe3_M2 + Ti_M2 + Mg_M2 + Fe2_M2;
table_cation.Total_M2 = Total_M2;


% Site occupancy indicators
% A site occupancy rate
XvA = vA;
table_cation.XvA = XvA;
% Na occupancy at A site
XNa_A = Na_A2;
table_cation.XNa_A = XNa_A;
% K occupancy at A site
XK_A = K_A2;
table_cation.XK_A = XK_A;
% Ca occupancy at M4 site
% Divide by 2 because the maximum capacity of M4 site is 2
XCa_M4 = Ca_M4 ./ 2;
table_cation.XCa_M4 = XCa_M4;
% Na occupancy at M4 site
XNa_M4 = Na_M4 ./ 2;
table_cation.XNa_M4 = XNa_M4;
% cm occupancy at M4 site
Xcm = cm4 ./ 2;
table_cation.Xcm = Xcm;
% Mg occupancy at M1 & M3 site
XMg_M1M3 = Mg_M1M3 ./ 3;
table_cation.XMg_M1M3 = XMg_M1M3;
% Fe2+ occupancy at M1 & M3 site
XFe2_M1M3 = Fe2_M1M3 ./ 3;
table_cation.XFe2_M1M3 = XFe2_M1M3;
% Mg occupancy at M2 site
XMg_M2 = Mg_M2 ./ 2;
table_cation.XMg_M2 = XMg_M2;
% Fe2+ occupancy at M2 site
XFe2_M2 = Fe2_M2 ./ 2;
table_cation.XFe2_M2 = XFe2_M2;
% Fe3+ occupancy at M2 site
XFe3_M2 = Fe3_M2 ./ 2;
table_cation.XFe3_M2 = XFe3_M2;
% Al occupancy at M2 site
XAl_M2 = Al_M2 ./ 2;
table_cation.XAl_M2 = XAl_M2;
% Ti occupancy at M2 site
XTi_M2 = Ti_M2 ./ 2;
table_cation.XTi_M2 = XTi_M2;
% Al occupancy at T1 site
XAl_T1 = Al_T1 ./ 4;
table_cation.XAl_T1 = XAl_T1;
% Si occupancy at T1 site
XSi_T1 = Si_T1 ./ 4;
table_cation.XSi_T1 = XSi_T1;


% Site Occupancy Balance Calculation 1
% edenite component
ed = XNa_A;
table_cation.ed = ed;
% K-edenite component
ked = XK_A;
table_cation.ked = ked;
% magnesio-riebeckite
mr = XFe3_M2;
table_cation.mr = mr;
% cation minor
cm_balance = Xcm;
table_cation.cm_balance = cm_balance;
% Ti at M2 site
ti = XTi_M2;
table_cation.ti = ti;
% glaucophane component
gl = XNa_M4 - XFe3_M2;
table_cation.gl = gl;
% tschermakite component
ts = XAl_M2 - XNa_M4 + XFe3_M2;
table_cation.ts = ts;
% ferro-tschermakite component
ftr = XFe2_M2;
table_cation.ftr = ftr;
% tremolite component
tr = XMg_M2 - XNa_A - XK_A - cm_balance;
table_cation.tr = tr;
% total
Total_balance = ed + ked + mr + cm + ti + gl + ts + ftr + tr;
table_cation.Total_balance = Total_balance;

% Site Occupancy Balance Calculation 2
% ed2: Total occupancy rate of A site (Na + K)
ed2 = ed + ked;
table_cation.ed2 = ed2;
% mr: Fe3+ occupancy rate of M2 site
mr2 = mr;
table_cation.mr2 = mr2;
% gl2: glaucophane component
gl2 = gl;
table_cation.gl2 = gl2;
% ts2: the total amount of "oxidative" or "tetravalent alternative" 
% components at the M2 site
ts2 = ti + ts;
table_cation.ts2 = ts2;
% ftr2: Fe²⁺-derived reductive tschermakite component
ftr2 = ftr;
table_cation.ftr2 = ftr2;
% tr2: reductive site occupancy of M2 + M4 sites
tr2 = cm_balance + tr;
table_cation.tr2 = tr2;
% total
Total_balance2 = ed2 + mr2 + gl2 + ts2 + ftr2 + tr2;
table_cation.Total_balance2 = Total_balance2;

% Site occupancy balance indicators
% ked#: ratio of K at A site
kednum = ked ./ (ed + ked);
table_cation.kednum = kednum;
% ti#: ratio of Ti in tschermakite component
tinum = ti ./ (ti + ts);
table_cation.tinum = tinum;
% cm/tr
cm_tr = cm ./ tr;
table_cation.cm_tr = cm_tr;

end

%% (13) calcGrt
function table_cation = calcGrt(table_cation)
% Syntax:
% function table_cation = calcGrt(table_cation)
%
% Inputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Outputs:
% ------
% - table_cation @table
%     cation dataset calculated from EPMA dataset
%
% Calculate various parameters of Garnet
% Store the calculation results in "table_cation" and output them

% RFe
% RFe = Fe3+/(Fe3++Fe2+)
RFe = table_cation.Fe3_cation_normalized./(table_cation.Fe3_cation_normalized + table_cation.Fe_cation_normalized);
table_cation.RFe = RFe;

% For each element, 
% normalize by assuming the number of cations in the garnet to be 8.
% Si
Si_Grt = table_cation.Si_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Si_Grt = Si_Grt;
% Ti
Ti_Grt = table_cation.Ti_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Ti_Grt = Ti_Grt;
% Al
Al_Grt = table_cation.Al_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Al_Grt = Al_Grt;
% Fe2
Fe2_Grt = table_cation.Fe_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Fe2_Grt = Fe2_Grt;
% Mn
Mn_Grt = table_cation.Mn_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Mn_Grt = Mn_Grt;
% Mg
Mg_Grt = table_cation.Mg_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Mg_Grt = Mg_Grt;
% Ca
Ca_Grt = table_cation.Ca_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Ca_Grt = Ca_Grt;
% Na
Na_Grt = table_cation.Na_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Na_Grt = Na_Grt;
% K
K_Grt = table_cation.K_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.K_Grt = K_Grt;
% Cr
Cr_Grt = table_cation.Cr_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Cr_Grt = Cr_Grt;
% Ni
Ni_Grt = table_cation.Ni_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Ni_Grt = Ni_Grt;
% V
V_Grt = table_cation.V_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.V_Grt = V_Grt;
% P
P_Grt = table_cation.P_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.P_Grt = P_Grt;
% Fe3
Fe3_Grt = table_cation.Fe3_cation_normalized * 8 ./ table_cation.Total_cation_normalized;
table_cation.Fe3_Grt = Fe3_Grt;

% Total
Total_Grt = sum([Si_Grt, Ti_Grt, Al_Grt, Fe2_Grt, Mn_Grt, Mg_Grt, Ca_Grt, Na_Grt, K_Grt, Cr_Grt, Ni_Grt, V_Grt, P_Grt, Fe3_Grt], ...
    2, 'omitnan');
table_cation.Total_Grt = Total_Grt;

% Calculate the single component ratios of garnet
% Total of four components: Fe2+, Mn, Mg, and Ca
component_total = Fe2_Grt + Mn_Grt + Mg_Grt + Ca_Grt;
% almandin component
alm = Fe2_Grt ./ component_total;
table_cation.Alm = alm;
% spessartine component
sps = Mn_Grt ./ component_total;
table_cation.Sps = sps;
% pyrope component
prp = Mg_Grt ./ component_total;
table_cation.Prp = prp;
% grossular component
grs = Ca_Grt ./ component_total;
table_cation.Grs = grs;

end

%% (14) columns2Mean
function [table_data_mean] = columns2Mean(table_data)
% Syntax:
% function [table_data_mean] = columns2Mean(table_data)
%
% Inputs:
% ------
% - table_data: @table
%     EPMA dataset to be calculated into average values [mol]
%
% Outputs:
% ------
% - table_data_mean @table
%     table containing a result of calculation
%
% Average the calculated cation data for each sample code.

% Initialize a new table
table_data_mean = table();

% Check the size of "table_data"
s = size(table_data);
% table_dataの高さ
h = s(1);

% Delete the text to the right of "_" in the Comment column
% This procedure may generate multiple data with the same sample code
for i = 1:h
    sample_code = split(table_data.Comment(i),'_');
    table_data.Comment(i) = sample_code(1);
end

% Extract sample code from the "Comment" column without duplication
sample_code_unique = unique(table_data.Comment);

% Repeat for all sample codes
for i = 1:length(sample_code_unique)
    target_sample = sample_code_unique(i);

    % Extract data with the same sample code
    % logical corresponding to "target_sample"
    id = strcmp(table_data.Comment, target_sample);
    target_rows = table_data(id, :);

    % Extract cells that contain numeric data
    numeric_columns = target_rows(:, vartype('numeric'));

    % Calculate the mean value for each column
    mean_values = varfun(@mean, numeric_columns);

    % Remove the prefix "mean_" and revert to the original column names
    mean_values.Properties.VariableNames = strrep(mean_values.Properties.VariableNames, 'mean_', '');

    % Add the calculation result to the bottom
    table_data_mean = [table_data_mean; mean_values]; %#ok<AGROW>
end

% Add a sample code column to the left of the output table
table_data_mean = addvars(table_data_mean, sample_code_unique(:, 1), 'Before', 1, 'NewVariableNames', "SampleCode");

end
