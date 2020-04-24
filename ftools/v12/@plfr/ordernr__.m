function order_nr = ordernr__(pLFR,varname)
%% ordernr__
%
%  File: ordernr__.m
%  Directory: 7_ftools/ftools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 29. (2019b)
%

%%

if isfield(pLFR.Preferred_Variable_Order,varname)
    order_nr = pLFR.Preferred_Variable_Order.(varname);
else
    order_nr = pLFR.Preferred_Variable_Order.other;
end

end