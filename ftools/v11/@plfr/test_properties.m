function pLFR = test_properties(pLFR)
%%

failed = false;

for prop = properties(pLFR).'
    if strcmp(pLFR.(prop{1}),'[EMPTY]')
        fprintf('Property %s is empty\n', prop{1});
        failed = true;
    end
end

if failed
    pcz_dispFunctionStackTrace('first',0,'last',0);
end

end