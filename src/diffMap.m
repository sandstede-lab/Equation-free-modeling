%% The function runDiffMap takes in data from traffic simulations,
% builds a diffusion map, compares the eigen directions given by the
% diffusion map eigenvectors, and plots relevent eigenvectors scaled by
% corresponding eigenvalues
function [vec,val,eps] = runDiffMap(data,k, weight)

% build the data from the simulation results
allTime = data; % data should be passed in as headwaysD

% calcuate the pairwise distances between data points
D = squareform(pdist(allTime'));
eps = weight*median(D(:)); % choose epsilon for the kernel based on the pairwise distances

[vec,val] = diffusionMap(eps,D,k);          % calculate the diffusion map

%{
calculate how unique each eigen direction is
r = zeros(k, 1);
r(1) = 1;
for j = 2:k
    r(j) = linearFit(vec,j);
end
%}


end

    % plotEps creates a log-log plot of the number of data points that are
    % less than epsilon vs. epsilon
    function plotEps(distances)
        % create the values of epsilon to test from 0 to maxEps by stepSize
        epsilon = logspace(-3,2,1000);
        % count the number of points that are less than each epsilon
        L = zeros(size(epsilon));
        for iEps = 1:length(epsilon)
            curEps = epsilon(iEps);
            L(iEps) = sum(sum(distances<curEps,1),2);
        end
        % plot the log-log plot of L(epsilon) vs epsilon
        figure;
        loglog(epsilon,L);
        xlabel('\epsilon','FontSize',24);
        ylabel('L(\epsilon)','FontSize',20);
    end
    
    
%$ diffusionMap returns the diffusion map embedding of a given data set
% epsilon - the parameter to use in the affinity matrix
% distMatrix - matrix of pairwise distances between data points
% k - the number of eigenvectors to return
% Returns:
% vec - the first k nontrivial eigenvectors
% val - the first k nontrivial eigenvalues
function [vec,val] = diffusionMap(epsilon,distMatrix,k)
    %find Markov matrix and its eigenvalues/eigenvectors
    A = basicKernel(distMatrix);
    M = markovify(A);
    [eigenvec,eigenval] = eigs(M,k+1);
    
    %sort eigenvalues and eigenvectors
    seval = diag(sort(diag(eigenval),'descend'));
    [~,ind] = sort(diag(eigenval),'descend');
    sevec = eigenvec(:,ind);
    seval = sparse(seval);
    eigsigns = sevec(1,:)./abs(sevec(1,:));
    sevec = sevec * diag(eigsigns);
    
    %return first k eigenvalues/vectors
    val = seval(2:end,2:end);
    vec = sevec(:,2:end);    

end

    % create a Markov matrix
    function M = markovify(af)
        M = zeros(size(af));
        rsum = sum(af,2);
        for r = 1:length(rsum)
            M(r,:) = af(r,:)./rsum(r);
        end
    end

    % kernel function
    function af = basicKernel(s)
        af = exp(-s.^2/epsilon^2);
    end
    
% Function linearFit compares the eigendirection given by an eigenvector to
% the directions given by the previous eigenvectors
% evecs - matrix containing all of the eigenvectors
% k - the eigenvector to study
% Returns:
% res - the linear fit of eigenvector k from the k-1 previous eigenvectors
function res = linearFit(evecs,k)
m = size(evecs,1);                      % number of data points
PHI = [ones(m,1) evecs(:,1:k-1)];       % k-1 previous eigvectors in columns with column of 1s in front
phi = evecs(:,k);                       % kth eigenvector
pdists = pdist(PHI(:,2:end));           % pairwise distances between eigenvectors in PHI
ereg = median(pdists)/3;                % epsilon given by the median of the pairwise distances
sqpdists = squareform(pdists);          % mxm version of the pairwise distances

% find the linear approximation
approx = zeros(size(phi));              % initalize an nx1 empty vector
for i = 1:m
    curdists = sqpdists(i,:);           % find the distances corresponding to i
    dists = [curdists(1:i-1) , curdists(i+1:end)];  % remove the ith element
    W = diag(kernel(dists));            % diagonalize the distances
    curPHI = [PHI(1:i-1,:) ; PHI(i+1:end,:)];   % find the ith entries of PHI with (i,i) removed
    curphi = [phi(1:i-1,:) ; phi(i+1:end,:)];   % find the ith entries of phi with (i,i) removed
    
    curPHI_W = curPHI'*W;               % compute PHI'*W
    approx(i) = PHI(i,:)*(pinv(curPHI_W*curPHI)*(curPHI_W*curphi)); % solve y = phi*(PHI'*W*PHI)^-1*(PHI'*W*phi)
end

res = norm(phi-approx)/norm(phi);       % compute the difference between the approximation and phi

end
