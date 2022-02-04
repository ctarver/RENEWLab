classdef PrecoderBypass < Precoder
    %PRECODERBYPASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = PrecoderBypass(varargin)
            % Parse the inputs.
            vars = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            addParameter(vars, 'name', 'Bypass', @(x) any(validatestring(x,{'PrecoderBypass'})));
            addParameter(vars, 'required_domain', 'bypass', @(x) any(validatestring(x,{'bypass'})));
            addParameter(vars, 'required_fs', 122.88e6, validScalarPosNum);
            addParameter(vars, 'index', 1, validScalarPosNum);
            addParameter(vars, 'n_ant', 8, validScalarPosNum);
            parse(vars, varargin{:});
            obj.save_inputs_to_obj(vars);
        end
        
        function out = subclass_use(obj, in)
            % PrecoderBypass Subclass_use
            % Passes output on all the number of streams we expect
            [n_users, n_sym, n_scs] = size(in);
            assert(n_users == 1, 'Bypass assumes 1 user data will go on all antennas. Change p.n_users to 1!');
            out = zeros(obj.n_ant, n_sym, n_scs);
            for i_ant = 1:obj.n_ant
                out(i_ant, :, :) =  in(1, :, :);
            end
        end
        
        function update(obj, channel)
            % update does nothing for a bypass.
        end
        
        function report(obj)
            
        end
    end
end
