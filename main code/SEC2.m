% Load the data
data = load("Sample_JointRV_2024.mat");
X = data.XY(1, :);
Y = data.XY(2, :);

%%
% Calculate histogram bin edges using Scott's rule
[~, bins_X] = histcounts(X, 120,'Normalization', 'probability');
[~, bins_Y] = histcounts(Y, 120,'Normalization', 'probability');

% Calculate 2D histogram
[H, xedges, yedges] = histcounts2(X, Y, bins_X, bins_Y, 'Normalization', 'probability');

% Create meshgrid for plotting (using bin centers)
xcenters = (xedges(1:end-1) + xedges(2:end)) / 2;
ycenters = (yedges(1:end-1) + yedges(2:end)) / 2;
[xpos, ypos] = meshgrid(xcenters, ycenters);

% Plot using surf
figure;
surf(xpos, ypos, H');  % Use H directly, no need to flatten or filter
xlabel('X');
ylabel('Y');
zlabel('Probability');
title('Joint Probability Distribution P(X, Y)');
view(3); % Set 3D view

%% Marginal PDF Calculation
% Check RV type for X
rvTypeX = RVType(X);

if strcmp(rvTypeX, 'Discrete')
    % Estimate marginal PDF of X using histogram
    unique_x = unique(X);
    if length(unique_x) > 1
        [countsX, centersX] = hist(X, unique_x); % Use unique values as centers
        marginal_pdf_X = countsX / length(X);

        figure;
        bar(centersX, marginal_pdf_X);
        xlabel('X');
        ylabel('Probability');
        title('Marginal Probability Distribution of X (Histogram)');
    else
        warning('Cannot calculate marginal PDF of X using histogram. Too few unique values.');
    end

elseif strcmp(rvTypeX, 'Continuous')
    % Estimate marginal PDF of X using ksDen
    n_x = length(X);
    bandwidth_x = std(X) * n_x^(-1/5); % Example bandwidth
    [marginal_pdf_X, centersX] = ksDen(X, bandwidth_x);

    figure;
    plot(centersX, marginal_pdf_X);
    xlabel('X');
    ylabel('Probability Density');
    title('Marginal Probability Density of X (KDE)');
end

% Check RV type for Y
rvTypeY = RVType(Y);

if strcmp(rvTypeY, 'Discrete')
    % Estimate marginal PDF of Y using histogram
    unique_y = unique(Y);
    if length(unique_y) > 1
        [countsY, centersY] = hist(Y, unique_y); % Use unique values as centers
        marginal_pdf_Y = countsY / length(Y);

        figure;
        stem(centersY, marginal_pdf_Y);
        xlabel('Y');
        ylabel('Probability');
        title('Marginal Probability Distribution of Y (Histogram)');
    else
        warning('Cannot calculate marginal PDF of Y using histogram. Too few unique values.');
    end

elseif strcmp(rvTypeY, 'Continuous')
    % Estimate marginal PDF of Y using ksDen
    n_y = length(Y);
    bandwidth_y = std(Y) * n_y^(-1/5); % Example bandwidth
    [marginal_pdf_Y, centersY] = ksDen(Y, bandwidth_y);

    figure;
    plot(centersY, marginal_pdf_Y);
    xlabel('Y');
    ylabel('Probability Density');
    title('Marginal Probability Density of Y (KDE)');
end


%%
function rvType = RVType(RV)
    uniqueValues = unique(RV);
    ratioUnique = length(uniqueValues) / length(RV);
    threshold = 0.05;  % Adjust threshold as needed
    if ratioUnique > threshold
        rvType = 'Continuous';
    else
        rvType = 'Discrete';
    end
    fprintf('Ratio of unique values for %s: %.4f\n', inputname(1), ratioUnique);
    fprintf('Determined RV type: %s\n', rvType);

end

function [pdf, x] = ksDen(data, bandwidth)
    min_x = min(data) - 3 * bandwidth;
    max_x = max(data) + 3 * bandwidth;
    x_points = ceil(sqrt(length(data)));  % Number of points for evaluation
    x = linspace(min_x, max_x, x_points);
    pdf = zeros(size(x));
    n = length(data);
    for i = 1:length(x)
        for j = 1:n
            pdf(i) = pdf(i) + gaussianKernel((x(i) - data(j)) / bandwidth);
        end
        pdf(i) = pdf(i) / (n * bandwidth);
    end
end

function K = gaussianKernel(x)
    K = (1 / sqrt(2 * pi)) * exp(-0.5 * x.^2);
end