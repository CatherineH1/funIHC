function [label,clusters,fd,mean_clusters_mat] = funIHC(x,Data,del,type)

if nargin<3
    del=0.01;
end 

if nargin<3
    type = 0;
end 

%  set up the basis with a knot at every year
nT       = length(x);  
nC       = size(Data,2);
rng      = [min(x),max(x)];
nbasis   = nT-5 + 4;
norder   = 6;
basisobj = create_bspline_basis(rng, nbasis, norder);

%  smooth the data by penalizing the second derivative
options = optimset('LargeScale',  'off',  ...
                   'Display',     'off',   ...
                   'GradObj',     'off',     ...
                   'Hessian',     'off',    ...               
                   'TolX',        del,...
                   'PlotFcns','optimplotfval');

Lfdobj  = 2;
GCV_fun = @(lambda) smooth_basis_gcv(lambda, x, Data, basisobj, Lfdobj);
[optlambda,fval,exitflag,output] = fminbnd(GCV_fun, 10.^-6, 10^6, options);


fdParv  = fdPar(basisobj, Lfdobj, optlambda);
[fd, df, ~, ~, SSE]  = smooth_basis(x, Data, fdParv);    

stderr = sqrt(sum(SSE)/(nC*(nT-df)));  

disp(['standard error of fit = ',num2str(stderr)])

xfine   = linspace(x(1),x(end),101)';    

if(type==0)
C       = eval_fd(xfine,fd,0);
elseif(type==1)
C       = eval_fd(xfine,fd,1);   
elseif(type==2)
C       = getcoef(fd);
end

%Choose Alpha
PCorR = reshape(1-pdist2(C',C','Spearman'),nC.*nC,1);
if(any(PCorR<0))
PCorR(PCorR<0)=[];
end
if(any(PCorR>0.99999))
PCorR(PCorR>0.99999)=[];
end
alpha = linspace(min(PCorR)+1e-6, max(PCorR)-1e-6,20);

for i=1:length(alpha)
label  = IHC(C',alpha(i));
EC(i)  = evalclusters(Data',label,'silhouette').CriterionValues;
NCL(i) = max(label);
if(NCL(i)>=(nC))
    break;
end
end

%Remove Nans (too low) and zeros (too high)
ind_nan=isnan(EC);
if(~isempty(ind_nan))
EC=EC(~ind_nan);
NCL=NCL(~ind_nan);
alpha=alpha(~ind_nan);
end
alpha(knee_pt(NCL):end)=[];

alpha = linspace(alpha(1), alpha(end),20);

for i=1:length(alpha)
label  = IHC(C',alpha(i));
L    = tabulate(label);
DL{i} = L(:,2);
EC(i)  = evalclusters(Data',label,'silhouette').CriterionValues;
NCL(i) = max(label);
end

tmp=tabulate(NCL);
indT=find(ismember(NCL,tmp(find(tmp(:,2)==max(tmp(:,2))),1)));
rind = flip(indT);
[~,ind]=max(EC(rind));

%Optimal Alpha
[label,~,~,mean_clusters_mat,clusters]=IHC(C',alpha(rind(ind)));

end

function [label,fidxcluster,rmclusters,mean_clusters_mat,clusters]=IHC(data,alpha)

n = size(data,1);

%Cluster the data
idxclustertmp = clusterdata(data,'criterion','distance','cutoff',1-alpha,'distance','Spearman','linkage','average');

%Put the id in each cluster
for j = 1:max(idxclustertmp)
idxcluster{1}{j}=find(idxclustertmp==j);
end

orig_idxcluster{1} = idxcluster{1};

%Obtain the clusters
clusters      = cell(length(idxcluster{1}),1);
rmclusters    = cell(length(idxcluster{1}),1);
for j = 1:length(idxcluster{1})
    clusters{j} = data(orig_idxcluster{1}{j},:);
end

mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false); 
mean_clusters_mat = cell2mat(mean_clusters_mat);

k=2;
con =1;

if(size(mean_clusters_mat,1)==1)
    fidxcluster = orig_idxcluster{k-1};
    label = idxclustertmp;
    return
end

while(con)
    
     clear('idxclustertmp');
     
     if(size(mean_clusters_mat,1)==1)
         k=k-1;
         break;
     end

     %Cluster the centres of the clusters
     idxclustertmp = clusterdata(mean_clusters_mat,'criterion','distance','cutoff',1-alpha,'distance','Spearman','linkage','average');
     %idxclustertmp

     %Put the id in each cluster
     for j = 1:max(idxclustertmp)
     idxcluster{k}{j}=find(idxclustertmp==j);
     end
    
     clusters      = cell(length(idxcluster{k}),1);
     for j = 1:length(idxcluster{k})
         orig_idxcluster{k}{j}  = vertcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
         clusters{j} = data(orig_idxcluster{k}{j},:);
     end
     
    %% Prune
    for j = 1:length(clusters)
        PCorR = 1-pdist2(clusters{j},mod_mean(clusters{j}),'Spearman');
        %[j,any(PCorR<(alpha))]
        if(any(PCorR<(alpha)))
            indrm                = find(PCorR<alpha);
            rmclusters{k}{j}     = clusters{j}(indrm,:); %Save removed curves
            clusters{j}(indrm,:) = []; %remove curves from cluster
            rmclustersidx{k}{j}  = orig_idxcluster{k}{j}(indrm); %save id of removed curves
            orig_idxcluster{k}{j}(indrm) = []; %remove curves from id
        else
            rmclusters{k}{j} =[];
        end
    end
    
    if(~isempty(vertcat(rmclusters{k}{:})))  
          %put the removed curves into their own cluster
          orig_idxcluster{k}  = horzcat(orig_idxcluster{k},num2cell(vertcat(rmclustersidx{k}{:}))');
          clusters = vertcat(clusters,num2cell(vertcat(rmclusters{k}{:}),2));
    end

    mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
    mean_clusters_mat = cell2mat(mean_clusters_mat);
    
    PCorR = 1-pdist2(clusters{j},mod_mean(clusters{j}),'Spearman');

    %% Stopping Conditions
    
    %1. The indices of the clusters are identical
    %2. The number of iterations is >100
    %3. The number of jumping genes is less than 5%.
    
    for l = 1:length(orig_idxcluster{k})
    current_ind(orig_idxcluster{k}{l}) = l;
    end

    for p = 1:length(orig_idxcluster{k-1})
    previous_ind(orig_idxcluster{k-1}{p}) = p;
    end
    
    %Adjusted rand index   
    ARI(k) = randindex(current_ind, previous_ind);
        
    if(ARI(k)>0.999)
        break;
    end

    if(k>3)
    
    m = mean(abs(diff(ARI((k-3):k))));
         
    if((k>100) || m<0.02)
         
         con1 = 1;
         while(con1)
             
             k=k+1;
             
             if(size(mean_clusters_mat,1)==1)  
                 k=k-1;
                 break;
             end

             idxclustertmp = clusterdata(mean_clusters_mat,'criterion','distance','cutoff',1-alpha,'distance','Spearman','linkage','average');

             for j = 1:max(idxclustertmp)
             idxcluster{k}{j}=find(idxclustertmp==j);
             end
             
             clusters      = cell(length(idxcluster{k}),1);
             for j = 1:length(idxcluster{k})
                 orig_idxcluster{k}{j}  = vertcat(orig_idxcluster{k-1}{idxcluster{k}{j}});
                 clusters{j} = data(orig_idxcluster{k}{j},:);
             end
             mean_clusters_mat = cellfun(@mod_mean,clusters,'UniformOutput', false);
             mean_clusters_mat = cell2mat(mean_clusters_mat);
             
             C = 1-pdist(mean_clusters_mat,'Spearman');
             
             con1 = max(max(C))>(alpha);
      
         end
         break;
    end
   
    end
     
     %display(['No of cycles: ',num2str(k)]);

     k=k+1;
    
end

fidxcluster = orig_idxcluster{k};

for l = 1:length(fidxcluster)
    label(fidxcluster{l}) = l;
end
label=label';

end


function gcv = smooth_basis_gcv(lambda, x, Data, basisobj, Lfdobj)

fdParv  = fdPar(basisobj, Lfdobj, lambda);
[~, ~, gcv] = smooth_basis(x, Data, fdParv);

gcv=mean(gcv);
end

function[y] = mod_mean(x)

if(size(x,1)>1) 
    y = mean(x);
else
    y = x ;
end

end

function [res_x, idx_of_result] = knee_pt(y,x,just_return)
%set internal operation flags
use_absolute_dev_p = true;  %ow quadratic

%deal with issuing or not not issuing errors
issue_errors_p = true;
if (nargin > 2 && ~isempty(just_return) && just_return)
    issue_errors_p = false;
end

%default answers
res_x = nan;
idx_of_result = nan;

%check...
if (isempty(y))
    if (issue_errors_p)
        error('knee_pt: y can not be an empty vector');
    end
    return;
end

%another check
if (sum(size(y)==1)~=1)
    if (issue_errors_p)
        error('knee_pt: y must be a vector');
    end
    
    return;
end

%make a vector
y = y(:);

%make or read x
if (nargin < 2 || isempty(x))
    x = (1:length(y))';
else
    x = x(:);
end

%more checking
if (ndims(x)~= ndims(y) || ~all(size(x) == size(y)))
    if (issue_errors_p)
        error('knee_pt: y and x must have the same dimensions');
    end
    
    return;
end

%and more checking
if (length(y) < 3)
    if (issue_errors_p)
        error('knee_pt: y must be at least 3 elements long');
    end
    return;
end

%make sure the x and y are sorted in increasing X-order
if (nargin > 1 && any(diff(x)<0))
    [~,idx]=sort(x);
    y = y(idx);
    x = x(idx);
else
    idx = 1:length(x);
end

%the code below "unwraps" the repeated regress(y,x) calls.  It's
%significantly faster than the former for longer y's
%
%figure out the m and b (in the y=mx+b sense) for the "left-of-knee"
sigma_xy = cumsum(x.*y);
sigma_x  = cumsum(x);
sigma_y  = cumsum(y);
sigma_xx = cumsum(x.*x);
n        = (1:length(y))';
det = n.*sigma_xx-sigma_x.*sigma_x;
mfwd = (n.*sigma_xy-sigma_x.*sigma_y)./det;
bfwd = -(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det;

%figure out the m and b (in the y=mx+b sense) for the "right-of-knee"
sigma_xy = cumsum(x(end:-1:1).*y(end:-1:1));
sigma_x  = cumsum(x(end:-1:1));
sigma_y  = cumsum(y(end:-1:1));
sigma_xx = cumsum(x(end:-1:1).*x(end:-1:1));
n        = (1:length(y))';
det = n.*sigma_xx-sigma_x.*sigma_x;
mbck = flipud((n.*sigma_xy-sigma_x.*sigma_y)./det);
bbck = flipud(-(sigma_x.*sigma_xy-sigma_xx.*sigma_y) ./det);

%figure out the sum of per-point errors for left- and right- of-knee fits
error_curve = nan(size(y));
for breakpt = 2:length(y-1)
    delsfwd = (mfwd(breakpt).*x(1:breakpt)+bfwd(breakpt))-y(1:breakpt);
    delsbck = (mbck(breakpt).*x(breakpt:end)+bbck(breakpt))-y(breakpt:end);
    %disp([sum(abs(delsfwd))/length(delsfwd), sum(abs(delsbck))/length(delsbck)])
    if (use_absolute_dev_p)
        % error_curve(breakpt) = sum(abs(delsfwd))/sqrt(length(delsfwd)) + sum(abs(delsbck))/sqrt(length(delsbck));
        error_curve(breakpt) = sum(abs(delsfwd))+ sum(abs(delsbck));
    else
        error_curve(breakpt) = sqrt(sum(delsfwd.*delsfwd)) + sqrt(sum(delsbck.*delsbck));
    end
end

%find location of the min of the error curve
[~,loc] = min(error_curve);
res_x = x(loc);
idx_of_result = idx(loc);
end


function [ARI, RI] = randindex(labels1, labels2)
    % The function calculates the Rand index (RI) and Adjusted Rand index (ARI)
    % between two label assignments: labels1 and labels2.
    
    N = numel(labels1);  % Get the number of elements in the label vector
    
    % Initialize the four quantities: TP (true positive), FN (false negative), FP (false positive), TN (true negative)
    TP = 0; FN = 0; FP = 0; TN = 0;
    
    % Calculate TP, FN, FP and TN
    for i = 1:N-1
        for j = i+1:N
            if (labels1(i) == labels1(j)) && (labels2(i) == labels2(j))  % TP: Both labels1 and labels2 have the same class for samples i and j
                TP = TP + 1;
            elseif (labels1(i) == labels1(j)) && (labels2(i) ~= labels2(j))  % FN: labels1 have the same class for samples i and j, but labels2 not
                FN = FN + 1;
            elseif (labels1(i) ~= labels1(j)) && (labels2(i) == labels2(j))  % FP: labels2 have the same class for samples i and j, but labels1 not
                FP = FP + 1;
            else   % TN: Both labels1 and labels2 have different classes for samples i and j
                TN = TN + 1;
            end
        end
    end
    
    % Calculate Rand Index (RI)
    RI = (TP + TN) / (TP + FP + FN + TN);
    
    % Calculate Adjusted Rand Index (ARI)
    try
        C = confusionmat(labels1, labels2);  % Use built-in function to generate confusion matrix
    catch
        C = myConfusionMat(labels1, labels2);  % If built-in function is not available, use the self-defined function
    end
    
    % Compute necessary quantities for ARI calculation
    sum_C = sum(C(:));
    sum_C2 = sum_C * (sum_C - 1);
    sum_rows = sum(C, 2);
    sum_rows2 = sum(sum_rows .* (sum_rows - 1));
    sum_cols = sum(C, 1);
    sum_cols2 = sum(sum_cols .* (sum_cols - 1));
    sum_Cij2 = sum(sum(C .* (C - 1)));
    
    % Compute ARI
    ARI = 2 * (sum_Cij2 - sum_rows2 * sum_cols2 / sum_C2) / ...
        ((sum_rows2 + sum_cols2) - 2 * sum_rows2 * sum_cols2 / sum_C2);
end

function C = myConfusionMat(g1, g2)
    % This function generates a confusion matrix if the built-in function confusionmat is not available.
    % It takes two input vectors, g1 and g2, which represent two different label assignments.
    
    groups = unique([g1;g2]);  % Get the unique groups in g1 and g2
    
    C = zeros(length(groups));  % Initialize the confusion matrix with zeros
    
    % Calculate the confusion matrix
    for i = 1:length(groups)
        for j = 1:length(groups)
            C(i,j) = sum(g1 == groups(i) & g2 == groups(j));  % Count the number of samples that belong to group i in g1 and group j in g2
        end
    end
end

