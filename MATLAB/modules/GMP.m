classdef GMP < Module
    
    properties
        P
        M
        L
        use_conj
        use_dc
        use_even
        coeffs
        lambda
    end
    
    methods
        function obj = GMP(varargin)
            vars = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validBool = @(x) islogical(x);
            addParameter(vars, 'name', 'GMP', @(x) any(validatestring(x,{'GMP'})));
            addParameter(vars, 'required_domain', 'time', @(x) any(validatestring(x,{'time'})));
            addParameter(vars, 'required_fs', 122.88e6, validScalarPosNum);
            addParameter(vars, 'index', 1, validScalarPosNum);
            addParameter(vars, 'P', 7, validScalarPosNum);
            addParameter(vars, 'M', 4, validScalarPosNum);
            addParameter(vars, 'L', 0, validScalarPosNum);
            addParameter(vars, 'use_conj', false, validBool);
            addParameter(vars, 'use_dc', false, validBool);
            addParameter(vars, 'use_even', false, validBool);
            addParameter(vars, 'lambda', 0, validScalarPosNum);
            parse(vars, varargin{:});
            obj.save_inputs_to_obj(vars);
            
            % Create some default coeffs.
            obj.run_setup;
        end
        
        function y = use(obj, x)
            x = reshape(x, [], 1);
            X = obj.setup_basis_matrix(x);
            y = X * obj.coeffs;
        end
        
        function learn(obj, x, y)
            %learn will set up the GMP matrix and learn the coeffs.
            X = obj.setup_basis_matrix(x);
            obj.coeffs = obj.ls_estimation(X, y);
        end
        
        function beta = ls_estimation(obj, X, y)
            %ls_estimation
            % Solves problems where we want to minimize the error between a
            % lienar model and some input/output data.
            %
            %     min || y - X*beta ||^2
            %
            % A small regularlizer, lambda, is included to improve the
            % conditioning of the matrix.
            %
            
            % Trim X and y to get rid of 0s in X.
            X = X(obj.M+obj.L:end-obj.L, :);
            y = y(obj.M+obj.L:end-obj.L);
            
            beta = (X' * X + obj.lambda * eye(size((X' * X)))) \ (X' * y);
        end
        
        function X = setup_basis_matrix(obj, x)
            %setup_basis_matrix. Setup the basis matrix for the LS learning of
            %the PA parameters or for broadcasting through the PA model.
            %
            % obj.setup_basis_matrix(x)
            %
            % Inputs:
            %   x - column vector of the PA input signal.
            % Output:
            %   X - matrix where each column is the signal, delayed version of
            %   a signal, signal after going through a nonlinearity, or both.
            
            number_of_basis_vectors = numel(obj.coeffs);
            X = zeros(length(x), number_of_basis_vectors);
            
            if obj.use_even
                step = 1;
            else
                step = 2;
            end
            
            % Main branch
            count = 1;
            for i = 1:step:obj.P
                branch = x .* abs(x).^(i-1);
                for j = 1:obj.M
                    delayed_version = zeros(size(branch));
                    delayed_version(j:end) = branch(1:end - j + 1);
                    X(:, count) = delayed_version;
                    count = count + 1;
                end
            end
            
            % Lag term
            for k = 3:step:obj.P  % Lag/Lead doesn't exist for k=1
                absolute_value_part_base = abs(x).^(k-1);
                for m = 1:obj.L
                    lagged_abs = [zeros(m,1); absolute_value_part_base(1:end-m)];
                    main_base = x .* lagged_abs;
                    for l = 1:obj.M
                        X(l:end, count) = main_base(1:(end-l+1));
                        count = count + 1;
                    end
                end
            end
            
            % Lead term
            for k = 3:step:obj.P  % Lag/Lead doesn't exist for k=1
                absolute_value_part_base = abs(x).^(k-1);
                for m = 1:obj.L
                    lead_abs = [absolute_value_part_base(1+m:end); zeros(m,1)];
                    main_base = x .* lead_abs;
                    for l = 1:obj.M
                        X(l:end, count) = main_base(1:(end-l+1));
                        count = count + 1;
                    end
                end
            end
            
            if obj.use_conj
                % Conjugate branch
                for i = 1:step:obj.P
                    branch = conj(x) .* abs(x).^(i-1);
                    for j = 1:obj.M
                        delayed_version = zeros(size(branch));
                        delayed_version(j:end) = branch(1:end - j + 1);
                        X(:, count) = delayed_version;
                        count = count + 1;
                    end
                end
            end
            
            % DC
            if obj.use_dc
                X(:, count) = 1;
            end
        end
    end
    
    methods (Access = protected)
        function run_setup(obj)
            n_coeffs = obj.calc_number_of_order_terms * obj.M + ...
                2*((obj.calc_number_of_order_terms-1) * obj.M * obj.L);
            
            if obj.use_conj
                n_coeffs = 2*n_coeffs;
            end
            
            if obj.use_dc
                n_coeffs = n_coeffs + 1;
            end
            
            % Start MP coeffs being completely linear (no effect)
            obj.coeffs = zeros(n_coeffs, 1);
            obj.coeffs(1) = 1;
        end
        
        function number_of_terms = calc_number_of_order_terms(obj, order, even)
            if nargin == 1
                order = obj.P;
                even = obj.use_even;
            end
            
            if even
                number_of_terms = order;
            else
                number_of_terms = (order + 1) / 2;
            end
        end
    end
end
