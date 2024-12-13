function DVC=GetDefaultValuesAndConstants(S)

DVC= PrestoDefaultValuesAndConstants;
[b,Su]=GetUserDefinedParameters(S);
if (b==1) 
    warning('Presto can not get user-defined parameters');
end
if (b==2)
    DVC=AssignForStructures(DVC,Su);
end 
% $Id: GetDefaultValuesAndConstants.m,v 1.4 2002/08/30 16:05:00 fseyfert Exp $ 
