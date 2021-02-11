%$ diffusionMap returns the diffusion map embedding of a given data set
% epsilon - the parameter to use in the affinity matrix
% distMatrix - matrix of pairwise distances between data points
% k - the number of eigenvectors to return
% Returns:
% vec - the first k nontrivial eigenvectors
% val - the first k nontrivial eigenvalues
function [vec,val] = diffMap(epsilon,distMatrix,k)
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

end