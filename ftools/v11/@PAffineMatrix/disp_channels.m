function disp_channels(N)
%% disp_channels
%  
%  File: disp_channels.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%
    
    args_str = strjoin(cellfun(@(v) {char(v)}, num2cell(N.vars)),',');

    for i = 1:numel(N.channels)

        Thetai = N.get_matrices(N.channels(i));
        
        if isa(N.Theta,'sdpvar')
            fprintf('Values of %s(%s) for channel %s:\n', N.name, args_str, char(N.channels(i)))
            display(Thetai);
        else
            fprintf('Values of %s(%s) for channel %s:\n', N.name, args_str, char(N.channels(i)))
            disp(Thetai);
        end
    end

end