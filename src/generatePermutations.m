function [all] = generatePermutations()

    all = zeros(5^5, 5);
    t = 1;

    for i = 1:5
        for j = 1:5
            for k = 1:5
                for l = 1:5
                    for m = 1:5
                        all(t,:) = [i, j, k, l, m];
                        t = t + 1;
                    end;
                end;
            end;
        end;
    end;

return;