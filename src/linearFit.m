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

    % kernel function 
    function dist = kernel(d)
        dist = exp(-d.^2/ereg^2);
    end
end