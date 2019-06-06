clear, clc;

addpath(genpath('/ihome/ctseng/liz86/MOG/Code/SLEP'));
ind = csvread("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Data/ind.csv");

%----------------------- Set optional items -----------------------
opts=[];

% Starting point
opts.init=2;        % starting from a zero point

% Termination 
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations

% regularization
opts.rFlag=1;       % use ratio

% Normalization
opts.nFlag=0;       % without normalization

% Group Property
opts.ind= ind; 

opts.G=[1:300, 1:300];

z_vec = 0.01:0.01:1;
num_simu = 100;
BEST_Z = zeros(1, num_simu);
MSE = zeros(1, num_simu);
AUC = zeros(1, num_simu);

for s = 1:num_simu
	%----------------------- load Data -----------------------
	true_beta = csvread(strcat("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Data/U5/true_beta_", num2str(s), ".csv"));                                     
	X_train = csvread(strcat("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Data/U5/X_train_", num2str(s), ".csv"));
	Y_train = csvread(strcat("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Data/U5/Y_train_", num2str(s), ".csv"));
	X_test = csvread(strcat("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Data/U5/X_test_", num2str(s), ".csv"));
	Y_test = csvread(strcat("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Data/U5/Y_test_", num2str(s), ".csv"));

	%----------------------- analysis -----------------------
	inner_fold_id = csvread(strcat("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Data/U5/innerFoldID_", num2str(s), ".csv"));
	MSE_inner = zeros(length(z_vec),10);
	BETA_path = zeros(length(z_vec), size(X_train, 2)); 

	for z_index = 1:length(z_vec)
		z = z_vec(z_index);
		% inner cross-validation
		for testid = 1:10
			inner_train_id = find(inner_fold_id ~= testid);
			inner_test_id = find(inner_fold_id == testid);
			X_inner_train = X_train(inner_train_id, :);
			X_inner_test = X_train(inner_test_id, :);
			Y_inner_train = Y_train(inner_train_id);
			Y_inner_test = Y_train(inner_test_id);

			[beta_inner, funVal2, ValueL2]= tree_LeastR(X_inner_train, Y_inner_train, z, opts);
			Y_hat = X_inner_test*beta_inner;
			MSE_inner(z_index,testid) = mean(((Y_hat-Y_inner_test).^2));
		end

		% use entire traing set to get solution path
		[beta, funVal2, ValueL2]= tree_LeastR(X_train, Y_train, z, opts);
		BETA_path(z_index, :) = beta; 
	end

	% select best tuning parameter
	ave_mse = mean(MSE_inner, 2);
	best_z_index = find(ave_mse == min(ave_mse));
	BEST_Z(s) = z_vec(best_z_index);

	% get MSE corresponding to best tuning parameter
	beta_best = BETA_path(best_z_index,:);
	Y_hat = X_test*transpose(beta_best);
	MSE(s) = mean(((Y_hat-Y_test).^2));

	% get AUC
	selection_prob = mean(BETA_path~=0, 1);
	[Xp,Yp,Tp,AUCp] = perfcurve(true_beta~=0, selection_prob, 1);
	AUC(s) = AUCp;
end

csvwrite("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Results/U5_SLEP_BEST_Z.csv",BEST_Z);
csvwrite("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Results/U5_SLEP_MSE.csv",MSE);
csvwrite("/ihome/ctseng/liz86/MOG/Exp_revision1/Simu3_nonOverlapping/Results/U5_SLEP_AUC.csv",AUC);


mean(MSE) %18.4679
std(MSE) %6.0547
mean(AUC) %0.8823
std(AUC)  %0.0395
