% test cases for section 1

% Note: this files is saved in where the current matlab session is.
% you can know the path by typing "pwd" in the command line in command
% window.

% Test Case 2: X2 ~ U(-5, 2)
X2 = -5 + (2 - (-5)) * rand(1, 1000000);
save('test_case_2.mat', 'X2');

% Test Case 3: X3 ~ N(3, 4)
X3 = 3 + 2 * randn(1, 1000000);
save('test_case_3.mat', 'X3');

% Test Case 4: X4 ~ Bin(5, 0.3)
X4 = binornd(5, 0.3, 1, 1000000);
save('test_case_4.mat', 'X4');

% Test Case 5: X5 ~ Poisson(10)
X5 = poissrnd(10, 1, 1000000);
save('test_case_5.mat', 'X5');