classdef Sim_Array < Array
    %SIM_ARRAY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_antennas
    end
    
    methods
        function obj = Sim_Array(varargin)
            % Parse the inputs.
            vars = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0); 
            
            addParameter(vars, 'name', 'Sim_Array', @(x) any(validatestring(x,{'Sim_Array'})));
            addParameter(vars, 'required_domain', 'time', @(x) any(validatestring(x,{'time', 'freq'})));
            addParameter(vars, 'required_fs', 122.88e6, validScalarPosNum);
            addParameter(vars, 'index', 1, validScalarPosNum);
            addParameter(vars, 'n_antennas', 1, validScalarPosNum);
            parse(vars, varargin{:});
            obj.save_inputs_to_obj(vars);
        end
        
        function out = subclass_tx(obj, in)
            out = in;
        end
        
        function out = subclass_rx(obj, in)
            out = in;
        end
        
        function subscribe_rx(obj, in)
            % nop
        end
        
        function subclass_measure_noise(obj)
            
        end
        
        function report(obj)
            
        end
    end
end

