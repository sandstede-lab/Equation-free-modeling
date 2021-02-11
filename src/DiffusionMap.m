classdef DiffusionMap
    properties
        data
        distMatrix
        evecs
        evals
        epsilon
        lsqOptions
        
    end
    
    methods
        function obj = DiffusionMap(input_data, num_eigvecs, weight)
            % initialize object
            obj.data = input_data;
            obj.distMatrix = squareform(pdist(obj.data'));
            obj.epsilon =  weight*median(obj.distMatrix(:));
            
            %find Markov matrix and its eigenvalues/eigenvectors
            A = exp(-obj.distMatrix.^2/obj.epsilon^2);
            M = zeros(size(A));
            rsum = sum(A,2);
            for r = 1:length(rsum)
                M(r,:) = A(r,:)./rsum(r);
            end
            [eigenvec,eigenval] = eigs(M,num_eigvecs+1);
            
            %sort eigenvalues and eigenvectors
            seval = diag(sort(diag(eigenval),'descend'));
            [~,ind] = sort(diag(eigenval),'descend');
            sevec = eigenvec(:,ind);
            seval = sparse(seval);
            eigsigns = sevec(1,:)./abs(sevec(1,:));
            sevec = sevec * diag(eigsigns);
            
            %return first k eigenvalues/vectors
            obj.evals = seval(2:end,2:end);
            obj.evecs = sevec(:,2:end);
            
            % initialize lsqnonlin options
            obj.lsqOptions = optimset('Display','Off', 'TolX', 1e-12);
            obj.lsqOptions.FunctionTolerance = 1e-8;
            obj.lsqOptions.StepTolerance = 1e-9;
            obj.lsqOptions.OptimalityTolerance = 1e-8;
        end
        
        % restriction operator with the diffusion map
        % newData is a micro variable
        % pnew is the restriction of newData
        function pnew = restrict(obj, newData)
            dist = pdist2(newData',obj.data')';
            w = exp(-dist.^2/obj.epsilon^2);
            k = (1./sum(w)).*w;
            pnew = (obj.evecs' * k)./diag(obj.evals);
        end
        
        % lifting operator with the diffusion map
        % newVal is a macro variable
        % v0, h
        function lifted = lift(obj, newVal)
            
            % find closest datapoints
            validPoints = 10;
            newDist = pdist2(newVal',obj.evecs);
            [~, index] = sort(newDist, 2, 'ascend');
            idxMin = index(1:validPoints);
            closestPoints = obj.data(:, idxMin); 
            
            % solve for best coefficients for linear combination
            liftedCoeffs = lsqnonlin(@(coeffs)liftingEquations(obj, coeffs, closestPoints, newVal, validPoints), ...
                1/validPoints*ones(validPoints, 1), zeros(validPoints, 1), ones(validPoints,1), obj.lsqOptions);
            lifted = closestPoints * liftedCoeffs; % create new lifted profile
        end
        
        % sets up and solves a system of equations to find the best linear combination
        function out = liftingEquations(obj, coeffs, closestPoints, newVal, validPoints)
            liftGuess = closestPoints * coeffs;
            restrictGuess = restrict(obj, liftGuess);
            out = [restrictGuess - newVal  ; sum(coeffs) - 1; zeros(validPoints - 1, 1)];
        end
    end
    
end