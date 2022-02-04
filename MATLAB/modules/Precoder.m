classdef Precoder < Module
    %PRECODER Superclass for all Precoders
    
    properties
        n_ant
    end
    
    methods
        function obj = Precoder()
            %PRECODER Construct an instance of this class
        end
        
        function precoded_data = use(obj, S_in)
            % Inputs:
            % S             SuperSignal
            %
            % Outputs:
            % precoded_data SuperSignal
            
            S = S_in.copy();
            S.match_this(obj.required_domain)
            
            %% Use the selected precoding subclass
            precoded_out = obj.subclass_use(S.data);
            
            %% Pack the returned data in a mSignal object
            if strcmp(obj.required_domain, 'bypass')
                new_domain = S.domain;
            else
                new_domain = obj.required_domain;
            end
            
            precoded_data = Signal(precoded_out, obj.n_ant, new_domain, ...
                S.fs, S.modulator, 'Pre Data');
        end
    end
    
    methods (Abstract)
        subclass_use(obj, S);
        update(obj, H);
    end
end
