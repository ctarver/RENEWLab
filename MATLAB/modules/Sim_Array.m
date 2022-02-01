classdef Sim_Array < Module
    %SIM_ARRAY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_antennas
    end
    
    methods
        function obj = Sim_Array(varargin)
            % Parse the inputs.
            vars = inputParser;
            valid_constellations = {'BPSK', 'QPSK', '16QAM','64QAM'};
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validBool = @(x) islogical(x);
            
            
            addParameter(vars, 'name', 'Sim_Array', @(x) any(validatestring(x,{'Sim_Array'})));
            addParameter(vars, 'required_domain', 'time', @(x) any(validatestring(x,{'time'})));
            addParameter(vars, 'required_fs', 122.88e6, validScalarPosNum);
            addParameter(vars, 'index', 1, validScalarPosNum);
            addParameter(vars, 'n_antennas', 1, validScalarPosNum);
            parse(vars, varargin{:});
            
            % Save inputs to obj
            fields = fieldnames(vars.Results);
            for i = 1:numel(fields)
                obj.(fields{i}) = vars.Results.(fields{i});
            end
        end
        
        function tx(obj)
            
        end
        
        function rx(obj)
            
        end
        function report(obj)
            
        end
    end
end

