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
        function pnew = restrict(obj, newData, i)
            allData = obj.data;
            vecs = obj.evecs;
            
            if nargin > 2
                % remove this point from the data set
                allData(:, i) = [];
                vecs(i,:) = [];
            end
            
            dist = pdist2(newData',allData')';
            w = exp(-dist.^2/obj.epsilon^2);
            k = (1./sum(w)).*w;
            pnew = (vecs' * k)./diag(obj.evals);
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
        
        %nystrom error comparison, calculated the following way:
        %   'correct' coordinates are given by the evec coord
        %   approximated coordinates are given by diff map restrict, when the value
        %   is removed from the data set
        %   better approximated coordinates are given by diff map restrict, when
        %   the value is still included in the data set
        function [percentError, restricted, restrictDiff] = testRestrict(obj, numRestrict)
            if nargin < 2
                numRestrict = length(obj.evecs);
            end
            
            restricted = zeros(numRestrict, length(obj.evals));
            restrictDiff = zeros(numRestrict, 1);
            percentError = zeros(numRestrict,1);
            
            for i = 1:numRestrict
                if mod(i, 500)==0
                    disp(i);
                end
                
                restricted(i, :) = restrict(obj,obj.data(:,i), i); % restrict this point
                restrictDiff(i) = norm((obj.evecs(i, :) - restricted(i, :)));
                percentError(i) = restrictDiff(i)/norm(obj.evecs(i, :));
            end
        end
        
        function[percentError, restricted, restrictDiff] = testLift(obj, numRestrict)
            if nargin < 2
                numRestrict = length(obj.evecs);
            end
            
            restricted = zeros(numRestrict, length(obj.evals));
            restrictDiff = zeros(numRestrict, 1);
            percentError = zeros(numRestrict,1);
            
            for i=1:numRestrict
                if mod(i, 500)==0
                    disp(i);
                end
   
                lifted = lift(obj, obj.evecs(i, :)');  % lift the profile
               % restrictDiff(i) = norm( lifted - obj.data(:,i));
               % percentError(i) = abs(lifte)-obj.evecs(i,:))/abs(obj.evecs(i,:));
                restricted(i, :) = restrict(obj, lifted); % restrict the lifted profile
                restrictDiff(i) = norm(restricted(i, :) - obj.evecs(i,:));    % compute the difference from the original embedding
                percentError(i) = abs(restricted(i, :) -obj.evecs(i,:))/abs(obj.evecs(i,:));
            end
        end
        
    end
    
end