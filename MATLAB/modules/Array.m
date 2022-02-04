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
        fb_paths
        rx_paths
        tx_paths
        fb_mode_enable
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
            
            % Extract the data to pass into the subclass tx function
            S_matrix = S.extract_data();
            
            % TODO. Put this in the extract data method so it's automatic.
            % But it needs to only be done for the array.
            if(strcmpi(S.mod_settings.name, 'OFDM'))
                i_clip = S.mod_settings.clip_index;
                S_matrix = S_matrix(:, i_clip+1: end-i_clip);
            end
            
            obj.S_last = S;
            obj.S_last.mod_settings.clip_index = 0; % There is no clip index since we've extracted the data already.
            
            % Call subclass_tx, wrap into mSignal, and return
            S_out = mSignal(obj.subclass_tx(S_matrix), obj.n_tx_paths, obj.required_domain, ...
                obj.required_fs, obj.S_last.mod_settings);
            
        end
        
        function S_out = rx(obj)
            % RX() Calls subclass receiver
            % Outputs:
            % S_out             mSignal
            
            out = obj.subclass_rx();
            
            %             if any(abs(real(out)) > 0.5) || any(abs(imag(out)) > 0.5)
            %                 warning('Array output is > 0.5. In your last bit! Possibille clipping! Be Careful!!!!!!');
            %             end
            
            % Format into SuperSignal and return streams based on mode
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
            
            obj.S_out_last = mSignal(out, size(out,1), obj.required_domain, ...
                obj.required_fs, obj.S_last.mod_settings);
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