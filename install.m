function install
% Generated with Toolbox Extender https://github.com/ETMC-Exponenta/ToolboxExtender
dev = ProbabilityAnalysisM4Dev;
dev.test('', false);
% Post-install commands
cd('..');
ext = ProbabilityAnalysisM4Extender;
ext.doc;
% Add your post-install commands below