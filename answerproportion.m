% Analysis script for MATLAB
% TODO: Proportion correct as a function of disc overlap in different conditions

% Read data
data = readtable('C:/Users/somme/Downloads/Coc1ben21.dat', 'Delimiter', '\t', 'HeaderLines', 0);

% Rename columns
data.Properties.VariableNames{'Var6'} = 'discs';
data.Properties.VariableNames{'Var9'} = 'condition';
data.Properties.VariableNames{'Var26'} = 'answers';

% Select relevant columns
data_reduced = data(:, {'discs', 'condition', 'answers'});

% Calculate proportion correct for each disc overlap and condition using splitapply
grouped_data = findgroups(data_reduced.condition, data_reduced.discs);
proportion_correct = splitapply(@mean, data_reduced.answers, grouped_data);

% Transpose the result to match the dimensions
proportion_correct = proportion_correct';

% Unique conditions and discs
unique_conditions = unique(data_reduced.condition);
unique_discs = unique(data_reduced.discs);

% Create a meshgrid of conditions and discs
[ConditionMesh, DiscMesh] = meshgrid(unique_conditions, unique_discs);

% Reshape data for plotting
proportion_correct_table = table();
proportion_correct_table.condition = ConditionMesh(:);
proportion_correct_table.discs = DiscMesh(:);
proportion_correct_table.prop_correct = proportion_correct';

% Plot the data using bar plot
figure;
bar(proportion_correct_table.condition, proportion_correct_table.prop_correct, 'FaceColor', [0.5 0.8 1], 'EdgeColor', 'k');
xlabel('Condition');
ylabel('Proportion Correct');
title('Proportion Correct by Condition');
xticks(unique_conditions);
grid on;