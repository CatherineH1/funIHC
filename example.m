%% Example of running funIHC
% The fda matlab library is required, please download from https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/

DataSet = load('U1505.mat'); %Data from Table 2, M=15 sigma = 0.05
Data = DataSet.Data; %The curves
id = DataSet.id; %The curves cluster assignment (ground truth)

x = 1:15; %Number of time points (M=15 or 200 in the simulations)
Data = transpose(Data);

funihc = funIHC(x,Data,0.01,0); %Run funIHC, type refers to the choice of distance metric: curves (0), first derivative (1) or coefficients (2)
randindex(funihc,id); %Compute ARI



