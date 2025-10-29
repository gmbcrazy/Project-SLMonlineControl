function FOVLevelCorrPlots(AveResTBL, Param, SaveP2)

    Param.xLim = [];
    Param.yLim = [];

    % Define all pairings and labels for iteration
    plotConfigs = {
        {'TargetSpeedR', 'Response_SpeedR', 'rNonTargetSpeed-Vs-rTargetSpeed', 'Corr (NonTargetRes.,SpeedR)', 'TargetSpeedR'},
        {'TargetSensoryR', 'Response_SpeedR', 'rNonTargetSpeed-Vs-rTargetSensory', 'Corr (NonTargetRes.,SpeedR)', 'TargetSensoryR'},
        {'TargetSpeedR', 'Response_SensoryR', 'rNonTargetSensory-Vs-rTargetSpeed', 'Corr (NonTargetRes.,SensoryR)', 'TargetSpeedR'},
        {'TargetSensoryR', 'Response_SensoryR', 'rNonTargetSensory-Vs-rTargetSensory', 'Corr (NonTargetRes.,SensoryR)', 'TargetSensoryR'},
        {'TargetSensoryR', 'Response', 'ResNonTarget-Vs-rTargetSensory', 'NonTargetRes.', 'TargetSensoryR'},
        {'TargetSpeedR', 'Response', 'ResNonTarget-Vs-rTargetSpeed', 'NonTargetRes.', 'TargetSpeedR'}
    };

    papersizePX = [0 0 17 8];

    % Loop over all 6 figures
    for i = 1:length(plotConfigs)
        cfg = plotConfigs{i};
        targetField = cfg{1};
        responseField = cfg{2};
        fileName = cfg{3};
        yLabelTxt = cfg{4};
        xLabelTxt = cfg{5};

        figure;
        [~, r, p, h] = LuPairRegressPlot_Group_Cov(AveResTBL.(targetField), AveResTBL.(responseField), AveResTBL.Speed, AveResTBL.Group, Param);

        ylabel(h(1), yLabelTxt);
        xlabel(h(1), xLabelTxt);
        ylabel(h(2), yLabelTxt);
        xlabel(h(2), xLabelTxt);
        text(h(2).XLim(1), h(2).YLim(2), 'Speed regressed', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

        % Save figure in multiple formats
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));

        print(gcf, [SaveP2 fileName '.svg'], '-dsvg', '-painters');
        print(gcf, [SaveP2 fileName '.tif'], '-dtiffn', '-painters');
        saveas(gcf, [SaveP2 fileName '.fig'], 'fig');
    end

    close all;
end