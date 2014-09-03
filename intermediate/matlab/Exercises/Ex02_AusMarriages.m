%% Import the data into a table.
M = readtable('AusMarriages.dat', 'Delimiter', '\t');

%% Insert an additional variable representing the numeric date.
M.dn = datenum(M.Date, 'dd/mm/yyyy');

%% Plot marriages over time.
figure
plot(M.dn, M.Marriages)
datetick('x', 'mmm yy', 'keepticks')
xlabel('Date')
ylabel('# of Marriages')
title('\bf Australian Marriages, 2007')

%% Smooth the data using conv.
n = 25;
v = ones(1, n)/n;
M.Msmth = conv(M.Marriages, v, 'same');

%% Overlay the smooth data.
hold on
plot( M.dn, M.Msmth, 'Color', 'r', 'LineWidth', 2)
hold off