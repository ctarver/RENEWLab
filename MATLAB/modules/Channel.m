classdef Channel < Module
    %CHANNEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        H   % Channel matrix
        n_users
        n_ant
    end
    
    methods
        function obj = Channel()
        end
        
        function ue_rx = use(obj, downlink_Signal_in)
            % Make copies of data
            downlink_mSignal = downlink_Signal_in.copy();
            downlink_mSignal.match_this(obj.required_domain, obj.required_fs);
            
            ue_rx = Signal(obj.subclass_use(downlink_mSignal.data), obj.n_users,...
                obj.required_domain, obj.required_fs, downlink_Signal_in.modulator, 'Channel Out');
        end
    end
end

