
data = load("Sample_JointRV_2024.mat");
X = data.XY(1, :);
Y = data.XY(2, :);

Z = (2.*X)-1;
W = 2-(3.*Y);

zType = RVType(Z);
wType = RVType(W);

numBins = ceil(sqrt(length(Z))); % Adjust as needed

figure;
if strcmpi(zType, 'Discrete')
[zCounts, zEdges] = histcounts(Z, numBins, 'Normalization', 'probability');
zXValues = (zEdges(1:end-1) + zEdges(2:end)) / 2;
bar(zXValues, zCounts);
else
zBandwidth = 0.1; % Adjust bandwidth as needed
[zPDF, zx] = ksDen(Z, zBandwidth);
plot(zx, zPDF);
end
title('Probability Density Function (PDF) of Z');
xlabel('Z Values');
ylabel('Probability Density');

figure;
if strcmpi(wType, 'Discrete')
[wCounts, wEdges] = histcounts(W, numBins, 'Normalization', 'probability');
wXValues = (wEdges(1:end-1) + wEdges(2:end)) / 2;
bar(wXValues, wCounts);
else
wBandwidth = 0.1;  % Adjust bandwidth as needed
[wPDF, wx] = ksDen(W, wBandwidth);
plot(wx, wPDF);
end
title('Probability Density Function (PDF) of W');
xlabel('W Values');
ylabel('Probability Density');

% Calculate histogram bin edges using Scott's rule
[~, bins_X] = histcounts(Z, 150,'Normalization', 'probability');
[~, bins_Y] = histcounts(W, 150,'Normalization', 'probability');

% Calculate 2D histogram
[H, xedges, yedges] = histcounts2(Z, W, bins_X, bins_Y, 'Normalization', 'probability');

% Create meshgrid for plotting (using bin centers)
xcenters = (xedges(1:end-1) + xedges(2:end)) / 2;
ycenters = (yedges(1:end-1) + yedges(2:end)) / 2;
[xpos, ypos] = meshgrid(xcenters, ycenters);

% Plot using surf
figure;
surf(xpos, ypos, H');  % Use H directly, no need to flatten or filter
xlabel('Z');
ylabel('W');
zlabel('Probability');
title('Joint Probability Distribution P(Z, W)');
view(3); % Set 3D view









%% Helper Functions
function rvType = RVType(RV)
uniqueValues = unique(RV);
ratioUnique = length(uniqueValues) / length(RV);
threshold = 0.05;
if ratioUnique > threshold
rvType = 'Continuous';
else
rvType = 'Discrete';
end
fprintf('Ratio of unique values: %.4f\n', ratioUnique);
fprintf('Determined RV type: %s\n', rvType);
end

function [PDF, x] = ksDen(data, bandwidth)
min_x = min(data) - 3 * bandwidth;
max_x = max(data) + 3 * bandwidth;
x = linspace(min_x, max_x, ceil(sqrt(length(data))));
PDF = zeros(size(x));
n = length(data);
for i = 1:length(x)
for j = 1:n
PDF(i) = PDF(i) + gaussianKernel((x(i) - data(j)) / bandwidth);
end
PDF(i) = PDF(i) / (n * bandwidth);
end
end

function K = gaussianKernel(x)
K = (1 / sqrt(2 * pi)) * exp(-0.5 * x.^2);
end