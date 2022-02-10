classdef Array < Module
    %ARRAY Superclass for all MIMO arrays.
    %     Properties
    %         name
    %         required_domain
    %         required_fs
    %         amp_scaling
    %         S_last
    %         n_fb_paths
    %         n_tx_paths
    %         n_rx_paths
    %
    %     Abstract
    %         subclass_tx(obj, X);
    %         X_HAT = subclass_rx(obj);
    %
    %     Public Methods
    %         S_out = rx(obj)
    %
    %     Static Methods
    %         obj = create(p)
    
    properties
        amp_scaling
        S_last            % Stores the last transmitted mSignal
        S_out_last        % Stores the last (full) received mSignal
        n_tx_paths
        n_rx_paths
        rx_paths
        tx_paths
        max_amp
    end
    
    methods (Abstract)
        X_HAT = subclass_tx(obj, X);
        X_HAT = subclass_rx(obj);
        X_HAT = subclass_measure_noise(obj);
    end
    
    methods
        function obj = Array()
        end
        
        function S_out = tx(obj, S_in)
            % TX() Calls subclass transmitter
            % Inputs:
            % S_in             mSignal
            
            % Ensure signals match block requirements
            S = S_in.copy();
            S.match_this(obj.required_domain, obj.required_fs);
            S.normalize_to_this_amp(obj.max_amp);
            
            % Extract the data to pass into the subclass tx function
            S_matrix = S.data;
            obj.S_last = S;
            
            fprintf('Transmitting %d samples\n', length(S.data));

            % Call subclass_tx, wrap into mSignal, and return
            S_out = Signal(obj.subclass_tx(S_matrix), S_in.n_streams, obj.required_domain, ...
                obj.required_fs, obj.S_last.modulator);
        end
        
        function S_out = rx(obj, in)
            % RX() Calls subclass receiver
            % Outputs:
            % S_out             mSignal
            
            out = obj.subclass_rx(in.data);
            S_out = Signal(out, obj.n_antennas, obj.required_domain, ...
                    obj.required_fs, in.modulator);
        end
        
        function S_out = measure_noise(obj)
            
            input('Disconnect RX ADC Ports.\n');
            out = obj.subclass_measure_noise();
            input('Reconnect RX ADC Ports.\n');
            
            % Format into SuperSignal and return
            if(obj.fb_mode_enable)
                S_out_rx = mSignal(out(obj.rx_paths,:), obj.n_rx_paths, obj.required_domain, ...
                    obj.required_fs, obj.S_last.mod_settings);
                S_out_fb = mSignal(out(obj.fb_paths,:), obj.n_fb_paths, obj.required_domain, ...
                    obj.required_fs, obj.S_last.mod_settings);
                S_out = [S_out_rx S_out_fb];
            else
                S_out = mSignal(out(obj.rx_paths,:), obj.n_rx_paths, obj.required_domain, ...
                    obj.required_fs, obj.S_last.mod_settings);
            end
        end
    end
end