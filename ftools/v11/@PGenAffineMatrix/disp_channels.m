function disp_channels(N)
%% disp_channels
%  
%  File: disp_channels.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%
    
    args_str = strjoin(cellfun(@(v) {char(v)}, num2cell(N.vars)),',');

    for i = 1:numel(N.channels)

        Thetai = N.get_matrixk(i);
        
        fprintf('Values of %s(%s) for channel #%d: %s\n', N.name, args_str, i, char(N.channels(i)))
        
        if isa(N.Theta,'sdpvar')
            display(Thetai);
        else
            disp(Thetai);
        end
    end

end