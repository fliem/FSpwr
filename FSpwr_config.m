clc
cfg = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. General

cfg.outputPath = '/path/to/outputDir';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Which measure would you like to analyze
% set cfg.measure to
% 1: cortical thickness
% 2: cortical surface area
% 3: cortical volume
% 4: subcortical volume

cfg.measure = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. Choose test
% set cfg.test to
% 'paired': paired sample t-test
% 'two-sample': independent (two-sample) sample t-test

cfg.test= 'two-sample';
%cfg.test= 'paired';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. What do you want to do
% set cfg.analysisType to
% 1: A-priori: calculate sample size that is required to detect an effect
%     given are:
%         * difference between groups
%         * alpha level
%         * required power
%
% 2: Post-hoc: calculate achieved power
%     given are:
%         * difference between groups
%         * alpha level
%         * sample size
%
% 3: Sensitivity: calculate required effect size
%     given are:
%         * alpha level
%         * power
%         * sample size

cfg.analysisType = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV. What are your parameters
% state your parameters (depending on what you chose at III.)
% set
% cfg.alpha (e.g., 0.05)
% cfg.tail (1: one-tailed, 2: two-tailed)
% cfg.power (is 1-beta, e.g., 0.8)
% cfg.N (sample size per group)
% cfg.difference (absolute differences between groups, e.g., in mm)
%

cfg.alpha = .05;
cfg.tail = 2;
cfg.power = .8;
cfg.difference = 0.25;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V. Run analysis with FSpwr(cfg);

FSpwr(cfg);


