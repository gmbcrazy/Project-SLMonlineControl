

function Result = NNMFScoresFromNeuroTrace(FOVres,NeuroTrace,SpeedTrace,NonTarget)
% COMPUTETESTSCORESFORFOV
% Use trained NNMF components stored in FOVres(iFOV) to score testing trials
% and assemble grouped datasets for plotting and summary tables.
%
% Outputs:
%   FOVres      - updated struct with InfoTBL/DataVec fields
%   ScoreGroup  - 3x2 struct array (group x whisk), each with fields:
%                 Score, ScoreNonTarget, Score_sub, Time
%   SpeedGroup  - 3x2 cell array of speed values per group/whisk condition
%   maxScore    - 1x2 max of each score column (for y-lims)
%   TBL_rows    - struct with fields: all, nontarget, sub (tables for concatenation)

    % SpeedTrace = BehTrace{1};




    % Prealloc containers
    DataVec = [];
    Score = [];
    ScoreNontarget = [];
    Score_sub = [];
    clear DataVecTBL
    tempScore = [];



    % Speed-regression correction (apply to first component only)
    pinvW = pinv(FOVres.W(:, FOVres.ComI));
    Score = (pinvW * NeuroTrace)';
    Score_speedreg=Score;
    % Score_speedreg(:,1) = Score(:,1) - [ones(size(Score,1),1) SpeedTrace] * FOVres.SpeedReg.B;
    Score_speedreg(:,1) = Score(:,1) - SpeedTrace * FOVres.SpeedReg.B(2);


    pinvW_sub = pinv(FOVres.W_NonTarget(:, FOVres.ComI_NonTarget)); %#ok<NASGU>
    Score_sub = (pinvW_sub * NeuroTrace(NonTarget,:))';
    Score_sub_speedreg=Score_sub;
    % Score_sub_speedreg(:,1) = Score_sub(:,1) - [ones(size(Score_sub,1),1) SpeedTrace] * FOVres.SpeedReg_NonTarget.B;
    Score_sub_speedreg(:,1) = Score_sub(:,1) - SpeedTrace * FOVres.SpeedReg_NonTarget.B(2);

    Result.Score=Score;
    Result.Score_speedreg=Score_speedreg;

    Result.Score_sub=Score_sub;
    Result.Score_sub_speedreg=Score_sub_speedreg;


end

