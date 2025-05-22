rng('default'); % For reproducibility

%% X2 ~ N(3, 4), Y2 ~ N(-5, 2) (Independent)
X = normrnd(3, 2, 1, 1000000);  % Generate X as a row vector
Y = normrnd(-5, sqrt(2), 1, 1000000); % Generate Y as a row vector
XY = [X; Y]; % Combine into a 2x1000000 matrix
save('case2.mat', 'XY');

%% TX3 ~ Gamma(2, 10), Y3 ~ Bin(4, 0.5) (Independent)
X = gamrnd(2, 10, 1, 1000000);  % Generate X as a row vector
Y = binornd(4, 0.5, 1, 1000000); % Generate Y as a row vector
XY = [X; Y]; % Combine into a 2x1000000 matrix
save('case3.mat', 'XY');

%% X4 ~ Exp(0.05), Y4 = 3X4 + 2 (Dependent)
X = exprnd(1/0.05, 1, 1000000);
Y = 3 * X + 2; 
XY = [X; Y];
save('case4.mat', 'XY');

%% TX5{-1, 1} (uniform) and Y5 = X5 + n, n ~ N(0, 0.5) (Dependent)
X = (rand(1, 1000000) < 0.5) * -1 + (rand(1, 1000000) >= 0.5) * 1; 
n = normrnd(0, sqrt(0.5), 1, 1000000); 
Y = X + n; 
XY = [X; Y]; 
save('case5.mat', 'XY');

disp('Test cases saved individually as separate .mat files.');