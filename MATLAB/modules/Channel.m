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
        
        function [enb_rx, ue_rx] = use(obj, downlink_mSignal_in, uplink_mSignal_in)
            % Make copies of data
            downlink_mSignal = downlink_mSignal_in.copy();
            uplink_mSignal = uplink_mSignal_in.copy();
            
            downlink_mSignal.match_this(obj.required_domain, obj.required_fs);
            uplink_mSignal.match_this(obj.required_domain, obj.required_fs);
            
            % Unpack the data.
            
            [enb_rx, ue_rx] = obj.subclass_use(X, Y);
        end
    end
end

