function [DeltaRatio, p] = PSTHtrace_PostPreTest(PSTHTrace, BaseI, ResI, TestMethod)

    Ref = PSTHTrace(:, BaseI);
    Res = PSTHTrace(:, ResI);

    % Baseline Normalization
    BaseNorm = min(Res(:));
    Ref = Ref - BaseNorm;
    Res = Res - BaseNorm;

    % Calculate Delta Ratio
    b = mean(Ref(:));
    r = mean(Res(:));
    DeltaRatio = (r - b) / b;

    % Perform Statistical Test
    switch lower(TestMethod)
        case 'ttest'
            error(['Test method ' TestMethod ' is not available.']);
        case 'ttest2'
            [H, p] = ttest2(Res(:), Ref(:));
        case 'ranksum'
            p = ranksum(Res(:), Ref(:));
        case 'signrank'
            error(['Test method ' TestMethod ' is not available.']);
        otherwise
            error(['Test method ' TestMethod ' is not available.']);
    end
end
