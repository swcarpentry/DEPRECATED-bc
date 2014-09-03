%% Import data from a single file containing heterogeneous data.
% We will use a table to hold the results.
medData = readtable('MedicalData.txt', 'Delimiter', '\t');
medData.Age = medData.Age/12; % Note: age is recorded in months.
medData.BPDiff = medData.BPSyst1-medData.BPDias1;
figure
scatter(medData.Age, medData.BPDiff)
xlabel('Age')
ylabel('Pressure')

%% Data from a file containing numeric data.
fileID = fopen('HeightWaistData.txt');
dataFormat = '%f %f'; % Two columns of numeric (double, floating point) data.
HWData = textscan(fileID, dataFormat, 'Headerlines', 1, 'Delimiter', '\t');
fclose(fileID);

% The data is imported into a 1x2 cell array (each column of data is stored
% in a separate cell).
figure; cellplot(HWData)
HW = cell2mat(HWData);

% Remove any observations containing NaNs.
Height = HW(:, 1);
Waist = HW(:, 2);
badObs = isnan(Height) | isnan(Waist);
cleanHWData = HW(~badObs, :);
disp(mean(cleanHWData))

%% What about reading data from multiple files?
fileList = cellstr(ls(['ArmsLegs', filesep, '*.txt']));

for k = numel(fileList):-1:1
    ArmsLegs(k).FileName = fileList{k};
    ArmsLegs(k).Data = dlmread(fileList{k}, '\t', 1, 0);
end

% Merge the data into one variable.
ALData = vertcat(ArmsLegs.Data);

females = strcmp(medData.Sex, 'F');
figure
scatter3(ALData(females, 1), ALData(females, 2), ALData(females, 3), 'k.')
hold on
scatter3(ALData(~females, 1), ALData(~females, 2), ALData(~females, 3), 'g.')
hold off
xlabel('Leg Length')
ylabel('Arm Length')
zlabel('Arm Circ')
legend('Females', 'Males')

%% Save data to a MAT-file.
save('S01_HealthData.mat')