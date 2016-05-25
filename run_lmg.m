contents=dir;

for i=1:size(contents,1)
        if isstrprop(contents(i).name(1), 'digit')&& contents(i).isdir==1
            cd(contents(i).name)
            lmg
            cd ..
        end
end

lmgAnalysis
    