coolness = round(rand(30, 1) * 100);
condition = repelem(1:5, 6);  % Create an array with each number repeated 6 times
condition = condition(randperm(length(condition)))';  % Randomize the order

unique_conditions = unique(condition);                % Get the unique condition values
avg_coolness = arrayfun(@(x) mean(coolness(condition == x)), unique_conditions);

% Plot the results
figure;
bar(unique_conditions, avg_coolness);                 % Bar plot of averages
xlabel('Condition');                                  % Label for x-axis
ylabel('Average Coolness');                           % Label for y-axis
title('Average Coolness per Condition');             % Title
grid on;    