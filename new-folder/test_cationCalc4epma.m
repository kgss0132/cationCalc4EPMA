% test_cationCalc4epma.m
% ------------------------------------------------------------------------
% prepare test file (input)
input_name = "testdata.xlsx";
input = readtable(input_name);
% ------------------------------------------------------------------------
% try cationCalc4epma
% please select "testdata.xlsx" from list dialog
% please select all minerals as target phases
% please select all 'Mass' data except for 'Total_Mass'
cationCalc4epma_eng_250925();

% ------------------------------------------------------------------------
% check (1): testdata_Selected.xlsx
% name of output file from cationCalc4epma
output_name1 = 'testdata_Selected.xlsx';
% name of expected file from cationCalc4epma
expected_name1 = 'expected_output_Selected.xlsx';
% sheet name of expected file
sheet_names1 = sheetnames(expected_name1);

% check all output files
for i = 1:numel(sheet_names1)
    % read output sheet
    output1 = readtable(output_name1,'Sheet',sheet_names1(i));
    % read expected sheet
    expected1 = readtable(expected_name1,'Sheet',sheet_names1(i));
    % compare two sheets ('Comment' raw)
    if ~isequaln(output1(:,:), expected1(:,:))
        error(['Test failed: values in table do not match expected output. Please see ' sheet_names1{i} ' sheet in ' output_name1]);
    end
end

% ------------------------------------------------------------------------
% check (2): testdata_Divided.xlsx
% name of output file from cationCalc4epma
output_name1 = 'testdata_Divided.xlsx';
% name of expected file from cationCalc4epma
expected_name1 = 'expected_output_Divided.xlsx';
% sheet name of expected file
sheet_names1 = sheetnames(expected_name1);

% check all output files
for i = 1:numel(sheet_names1)
    % read output sheet
    output1 = readtable(output_name1,'Sheet',sheet_names1(i));
    % read expected sheet
    expected1 = readtable(expected_name1,'Sheet',sheet_names1(i));
    % compare two sheets ('Comment' raw)
    if ~isequaln(output1(:,:), expected1(:,:))
        error(['Test failed: values in table do not match expected output. Please see ' sheet_names1{i} ' sheet in ' output_name1]);
    end
end

% ------------------------------------------------------------------------
% check (3): testdata_Results.xlsx
% name of output file from cationCalc4epma
output_name1 = 'testdata_Results.xlsx';
% name of expected file from cationCalc4epma
expected_name1 = 'expected_output_Results.xlsx';
% sheet name of expected file
sheet_names1 = sheetnames(expected_name1);

% check all output files
for i = 1:numel(sheet_names1)
    % read output sheet
    output1 = readtable(output_name1,'Sheet',sheet_names1(i));
    % read expected sheet
    expected1 = readtable(expected_name1,'Sheet',sheet_names1(i));
    % compare two sheets ('Comment' raw)
    if ~isequaln(output1(:,:), expected1(:,:))
        error(['Test failed: values in table do not match expected output. Please see ' sheet_names1{i} ' sheet in ' output_name1]);
    end
end

% ------------------------------------------------------------------------
% check (4): testdata_Mean.xlsx
% name of output file from cationCalc4epma
output_name1 = 'testdata_Mean.xlsx';
% name of expected file from cationCalc4epma
expected_name1 = 'expected_output_Mean.xlsx';
% sheet name of expected file
sheet_names1 = sheetnames(expected_name1);

% check all output files
for i = 1:numel(sheet_names1)
    % read output sheet
    output1 = readtable(output_name1,'Sheet',sheet_names1(i));
    % read expected sheet
    expected1 = readtable(expected_name1,'Sheet',sheet_names1(i));
    % compare two sheets ('Comment' raw)
    if ~isequaln(output1(:,:), expected1(:,:))
        error(['Test failed: values in table do not match expected output. Please see ' sheet_names1{i} ' sheet in ' output_name1]);
    end
end

% ------------------------------------------------------------------------
% report the completement of test
disp('================================================================')
disp('All output files are consistent with expected files!')
disp('================================================================')