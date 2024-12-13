function [b,P]=GetUserDefinedParameters(S)
b=0;
P=0;
if (isfield(S,'params'))
    b=1;
    try
        if (isstr(S.params))
            eval(S.params);
            P=DVC;
        else
            P=feval(S.params);
        end
            b=2;
    catch
    end
end
     
% $Id: GetUserDefinedParameters.m,v 1.3 2002/08/30 16:05:00 fseyfert Exp $ 
