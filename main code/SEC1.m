%1: unknown RV 
%2: X2 ~ U(-5, 2) 
%3: X3 ~ N(3, 4) 
%4: X4 ~ Bin(5, 0.3) 
%5: X5 ~ Poisson(10)

data = load("test_case_3.mat");
RV = data.X3;
rvType = RVType(RV);
n = length(RV);
numBins = ceil(sqrt(n));
tmax = 0.5;

if strcmpi(rvType, 'Discrete')
    % PDF/PMF and CDF Estimation
    [counts, edges] = histcounts(RV, numBins, 'Normalization','probability');
    xValues = (edges(1:end-1) + edges(2:end)) / 2; 
    
    % Plot PDF
    figure;
    bar(xValues, counts);
    title('Probability Density Function (PDF)');
    xlabel('Data Values');
    ylabel('Probability Density');
    
    % Calculate and plot CDF
    cdfValues = cumsum(counts) * (edges(2) - edges(1));
    figure;
    stairs(edges(1:end-1), cdfValues);
    title('Cumulative Distribution Function (CDF)');
    xlabel('x');
    ylabel('F(x)');
    grid on;
        

    % MGF, First Moment, and Second Moment Calculation
    t = linspace(0, tmax, 1000); 
    mgf = exp(t' * xValues) * counts'; 
    firstMoment = exp(t' * xValues) * (xValues .* counts)';
    secondMoment = exp(t' * xValues) * (xValues.^2 .* counts)';
    
    % Calculate first and second moment at t = 0
    firstMoment_at_0 = sum(counts .* xValues);  % E[X]
    secondMoment_at_0 = sum(counts .* (xValues.^2));  % E[X^2]
    
    % Plot MGF, First Moment, and Second Moment in subplots
    figure;
    subplot(3, 1, 1);
    plot(t, mgf, 'k');
    title('Moment Generating Function (MGF)');
    xlabel('t');
    ylabel('MGF');

    subplot(3, 1, 2);
    plot(t, firstMoment, 'r');
    title('First Moment');
    xlabel('t');
    ylabel('First Moment');
    
    subplot(3, 1, 3);
    plot(t, secondMoment, 'b');
    title('Second Moment');
    xlabel('t');
    ylabel('Second Moment');

    meanRV = mean(RV);
    varianceRV = var(RV);
    thirdMoment = meanRV^3;  % Calculate 3rd central moment
    
    % Display the results
    fprintf('Mean: %.4f\n', meanRV);
    fprintf('Variance: %.4f\n', varianceRV);
    fprintf('Third central Moment: %.4f\n', thirdMoment);
    fprintf('First Moment at t=0: %.4f\n', firstMoment_at_0);
    fprintf('Second Moment at t=0: %.4f\n', secondMoment_at_0);
    fprintf('Sum of PMF counts: %.4f\n', sum(counts));


elseif strcmpi(rvType, 'Continuous')
    bandwidth = 0.1;
    [PDF, x] = ksDen(RV, bandwidth);
    
    % Plot the PDF
    figure;
    plot(x, PDF);
    title('Estimated PDF using my\_ksdensity');
    xlabel('Value');
    ylabel('Density');
    
    % Calculate CDF using trapezoidal integration
    CDF = cumtrapz(x, PDF);
    
    % Plot CDF
    figure;
    plot(x, CDF);
    title('Estimated CDF');
    xlabel('Value');
    ylabel('Cumulative Probability');
    
    % Calculate MGF, First Moment, and Second Moment
    t = linspace(0, tmax, 1000);
    mgf = zeros(size(t));
    firstMoment = zeros(size(t));
    secondMoment = zeros(size(t));
    
    for i = 1:length(t)
        mgf(i) = trapz(x, exp(t(i) * x) .* PDF);
        firstMoment(i) = trapz(x, x .* exp(t(i) * x) .* PDF);
        secondMoment(i) = trapz(x, (x.^2) .* exp(t(i) * x) .* PDF);
    end
    
    % Calculate first and second moment at t = 0
    firstMoment_at_0 = trapz(x, x .* PDF);
    secondMoment_at_0 = trapz(x, (x.^2) .* PDF);
    
    % Plot MGF, First Moment, and Second Moment
    figure;
    subplot(3, 1, 1);
    plot(t, mgf, 'k');
    title('Moment Generating Function (MGF)');
    xlabel('t');
    ylabel('MGF');

    subplot(3, 1, 2);
    plot(t, firstMoment, 'r');
    title('First Moment');
    xlabel('t');
    ylabel('First Moment');
    
    subplot(3, 1, 3);
    plot(t, secondMoment, 'b');
    title('Second Moment');
    xlabel('t');
    ylabel('Second Moment');

    meanRV = mean(RV);
    varianceRV = var(RV);
    thirdMoment = mean((RV-meanRV).^3);
    
    % Display the results
    fprintf('Mean: %.4f\n', meanRV);
    fprintf('Variance: %.4f\n', varianceRV);
    fprintf('Third central Moment: %.4f\n', thirdMoment);
    fprintf('First Moment at t=0: %.4f\n', firstMoment_at_0);
    fprintf('Second Moment at t=0: %.4f\n', secondMoment_at_0);
end

%% Helper Functions Section
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
