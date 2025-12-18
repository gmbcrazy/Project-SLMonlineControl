
clear all
DataPath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step3\awakeRefSpon\GroupSLM20-Nov-2025\';
load([DataPath 'FOVoutputs.mat'],'SaveFunDate')
csvPath = [SaveFunDate '\LME_FixedEffectsSummary_AllWindows_Groups_SCALED.csv'];

ProcessPar.GroupColor=[255 51 153;91 20 212;121 247 111]/255;
pThresh = 0.05;
NDataName={'deltaF','spks'};

spec.DataType  = 'deltaF';
spec.Window    = 'Win3';
spec.Condition = 'Sensory';
spec.Subset    = 'Sensory_All';
spec.Model     = 'Extended_scaled';
spec.Model     = 'Base_scaled';
spec.Colormap  = slanCM('seismic',64);
spec.ResponseVar = 'Response';    % or any of the 4 responses


% Only include these factors and their pairwise interactions
spec.Term = {'Speed_z','SpeedR_z', 'TargetSpeedR_z','TargetCellN','SensoryR_z', 'TargetSensoryR_z'};
% spec.Term = {};
spec.TargetColor = ProcessPar.GroupColor(1,:);
spec.NodeColor =[0 0 0];
% Optional fixed color limits
spec.Clim = [-0.002 0.002];
pThresh = 0.05;
plotCol = 'Estimate';


% figure;
% [Gout, theta] = LME_CircosFromCSV(csvPath, spec, plotCol, pThresh);
% 


PlotColSubSet={'NonSensory','Sensory'};
PlotColSubSetLabel={'SLM','SLM + Sensory'};

PlotColSubCondition={'All','G1','G2','G3'};
PlotColSubConditionLabel={'All','L','S','N'};

FigureWin={'1','2','3'};

SaveLME=[SaveFunDate 'LMESummary\'];
mkdir(SaveLME);

plotValue={'Estimate','tValue'};
pValuesVec=[0.05 0.01 0.001];
pValuesLabel={'P05','P01','P001'};

ModelVec = {'Base_scaled','Extended_scaled'};


pThresh = 0.05;

 modelStrVec={'Response ~ (SpeedR_z + SensoryR_z + TargetCellN_z) * (TargetSpeedR_z + TargetSensoryR_z) + (SpeedR_z +TargetSpeedR_z)* Speed_z + (1 | Session)',...
     'Response ~ (SpeedR_z + SensoryR_z + TargetCellN_z) * (TargetSpeedR_z + TargetSensoryR_z) + (SpeedR_z +TargetSpeedR_z)* Speed_z + AveDist_z + MinDist_z + (1 | Session)'};


for iPth=1:length(pValuesVec)
    pThresh = pValuesVec(iPth);
    SaveLME_Sub=[SaveLME pValuesLabel{iPth} '\'];
    mkdir(SaveLME_Sub);

    for iModel=1:length(ModelVec)
        spec.Model     = ModelVec{iModel};

        for iPlotCol=1:length(plotValue)
            plotCol=plotValue{iPlotCol};
            if strmatch(plotCol,'Estimate')
               plotLabel='Coeff';
               spec.Clim=[-0.003 0.003];
            else
               plotLabel=plotCol;
               spec.Clim=[-10 10];
            end
        
            for iData=1:length(NDataName)
                spec.DataType  = NDataName{iData};
                for iWin=1:length(FigureWin)
                    spec.Window    = ['Win' FigureWin{iWin}];
                    figure;
                    for iRow=1:length(PlotColSubSet)
                        for jRow=1:length(PlotColSubCondition)
                            if jRow>1
                            spec.TargetColor = ProcessPar.GroupColor(jRow-1,:);
                            else
                            spec.TargetColor = spec.NodeColor;
                            end
                    
                            
                    
                    
                            spec.Condition = PlotColSubSet{iRow};
                            spec.Subset    = [PlotColSubSet{iRow} '_' PlotColSubCondition{jRow}];
                            t3{iRow,jRow}=subplotLU(length(PlotColSubSet),length(PlotColSubCondition),iRow,jRow,[0.05 0.01 0.1 0.1 0.05 0.05]);

                            if strcmp(spec.Model,'Base_scaled')
                               modelStr =  modelStrVec{1};
                            elseif strcmp(spec.Model,'Extended_scaled')
                               modelStr =  modelStrVec{2};
                            else
                               modelStr ='';
                            end

                            [Gout, theta] = LME_CircosFromCSV(csvPath, spec, plotCol, pThresh);
                            colormap(spec.Colormap);
                            set(gca,'clim',spec.Clim );
            
                            if jRow==1
                               ylabel(PlotColSubSetLabel{iRow},'Color',[0 0 0]);
                            end
                            if iRow==2
                               xlabel(PlotColSubConditionLabel{jRow},'Color',spec.TargetColor);
                            end
                    
                        end
                    end
            
                    bar2=colorbar(t3{1,1});
                    bar2.Location='northoutside';
                    bar2.Position=[t3{1,1}.Position(1)+0.05,t3{1,1}.Position(2)+t3{1,1}.Position(4)+0.01, 0.1,0.02];
                    bar2.Label.String=plotLabel;
                    bar2.Ticks=union(spec.Clim,0);
            

                    subplot('Position',[0.2 0.92 0.4 0.08]);
                    text(1,1,['LME: ' modelStr],'HorizontalAlignment','left','VerticalAlignment','top','FontSize',6);
                    hold on;
                    plot(8,0,'ko','MarkerSize',10,'MarkerFaceColor','k');
                    text(8.5,0,'Significant covariate','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
                    plot([14 16],[0 0],'k-',LineWidth=4);
                    text(16.5,0,'Significant interaction','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
                    
                    set(gca,'xlim',[0 20],'ylim',[0 1.5]);
                    axis off


            
                    papersizePX=[0 0 length(PlotColSubCondition)*6+2 length(PlotColSubSet)*6+2];
                    set(gcf, 'PaperUnits', 'centimeters');
                    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                    print(gcf, [SaveLME_Sub spec.DataType '_' spec.Model '_' plotLabel '_' spec.Window  '.svg'], '-dsvg', '-painters');
                    print(gcf, [SaveLME_Sub spec.DataType '_' spec.Model '_' plotLabel '_' spec.Window  '.tif'], '-dtiffn', '-painters');
                    close all
                end
            end
        end
    end
end





for iPth=1:length(pValuesVec)
    pThresh = pValuesVec(iPth);
    SaveLME_Sub=[SaveLME pValuesLabel{iPth} '\'];
    mkdir(SaveLME_Sub);

    for iModel=1:length(ModelVec)
        spec.Model     = ModelVec{iModel};

        for iPlotCol=1:length(plotValue)
            plotCol=plotValue{iPlotCol};
            if strmatch(plotCol,'Estimate')
               plotLabel='Coeff';
               spec.Clim=[-0.003 0.003];
            else
               plotLabel=plotCol;
               spec.Clim=[-10 10];
            end
        
            for iData=1:length(NDataName)
                spec.DataType  = NDataName{iData};
                for iWin=1:length(FigureWin)
                    spec.Window    = ['Win' FigureWin{iWin}];
                    figure;
                    for iRow=1:length(PlotColSubSet)

                            spec.TargetColor = spec.NodeColor;

                            spec.Condition = PlotColSubSet{iRow};
                            spec.Subset    = [PlotColSubSet{iRow} '_' PlotColSubCondition{1}];
                            t3{iRow}=subplotLU(1,length(PlotColSubSet),1,iRow,[0.05 0.01 0.15 0.1 0.05 0.05]);

                            if strcmp(spec.Model,'Base_scaled')
                               modelStr =  modelStrVec{1};
                            elseif strcmp(spec.Model,'Extended_scaled')
                               modelStr =  modelStrVec{2};
                            else
                               modelStr ='';
                            end

                            [Gout, theta] = LME_CircosFromCSV(csvPath, spec, plotCol, pThresh);
                            colormap(spec.Colormap);
                            set(gca,'clim',spec.Clim );
            

                            xlabel(PlotColSubSetLabel{iRow},'Color',[0 0 0]);


                    

                    end
            
                    bar2=colorbar(t3{1});
                    bar2.Location='northoutside';
                    bar2.Position=[t3{1,1}.Position(1)+0.05,t3{1,1}.Position(2)+t3{1,1}.Position(4)+0.01, 0.12,0.02];
                    bar2.Label.String=plotLabel;
                    bar2.Ticks=union(spec.Clim,0);
            

                    subplot('Position',[0.2 0.88 0.4 0.08]);
                    text(2,1,['LME: ' modelStr],'HorizontalAlignment','left','VerticalAlignment','top','FontSize',6);
                    hold on;
                    plot(8,0,'ko','MarkerSize',10,'MarkerFaceColor','k');
                    text(8.5,0,'Significant covariate','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
                    plot([16 18],[0 0],'k-',LineWidth=4);
                    text(18.5,0,'Significant interaction','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
                    
                    set(gca,'xlim',[0 20],'ylim',[0 1.5]);
                    axis off


            
                    papersizePX=[0 0 length(PlotColSubSet)*8+4 8+4];
                    set(gcf, 'PaperUnits', 'centimeters');
                    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                    print(gcf, [SaveLME_Sub 'AllGroup_' spec.DataType '_' spec.Model '_' plotLabel '_' spec.Window  '.svg'], '-dsvg', '-painters');
                    print(gcf, [SaveLME_Sub 'AllGroup_' spec.DataType '_' spec.Model '_' plotLabel '_' spec.Window  '.tif'], '-dtiffn', '-painters');
                    close all
                end
            end
        end
    end
end




% modelStr = LME_ModelFromCSV(csvPath, spec)