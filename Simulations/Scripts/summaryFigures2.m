% Run from within directory with spreadsheets giving FBA results for 2D matrix
% Takes in all the data spreadsheets, extracts specific parameters of
% interest and then generates summary plot

% Code for 'Whole-cell energy modeling reveals quantitative changes of
% predicted energy flows in RAS mutant cancer cell lines'

% .. Author: - Phil Luthert & Christina Kiel 11/5/22

% FUNCTIONS CALLED
    % plot4Lines()
    % plot4LinesIndependent()
    % matrixCorr()
    % figureSummary()
    
% DATA REQUIRED
    % dataRaw and dataKcat
    % "ATPGTPaseSummary.xlsx"
    
% [1] PARAMETERS FROM FBA ANALYSIS RE-USED HERE

load('dataKcat');
load('dataRaw');
dataPC = dataKcat;

% lose the kcat column
dataPC = removevars(dataPC,{'kcat'});

% now select subset, default kcat corrected values, columns 9-12
dataPC = dataPC(:,9:12);
S = dataPC.Properties.VariableNames(1:4); 

Squeeze = [0.0001 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.25 0.5 -1];
Multiplier = [0.003125 0.00625 0.0125 0.025 0.05 0.1 0.2 0.4 0.8 1.6 3.2];

% [2] SOME CONSTANTS / HOUSEKEEPING

disp('running summaryFigures');
clear('paramSummary');  % helpful if re-running
clear('listing');       % helpful if re-running
column = -1;
fileNames = S;
reshapeRow = size(Squeeze,2);

% Glucose uptake rates in pmol/ng protein/day from Table S1 Figure2F data
Qglc(1) = 43.96;      
Qglc(2) = 60.92;
Qglc(3) = 42.41;
Qglc(4) = 57.16;

summaryRows = size(Squeeze,2);
summaryCols = size(Multiplier,2);
sampleNo = size(S,2);

% [3] ALLOCATE SPACE etc

figureSummaryEXglcStore = zeros(summaryRows, summaryCols, sampleNo);
figureSummaryATPyieldStore = zeros(summaryRows, summaryCols, sampleNo);
figureSummaryDMatpStore = zeros(summaryRows, summaryCols, sampleNo);
figureSummaryPFKStore = zeros(summaryRows, summaryCols, sampleNo);
figureSummaryATPS4miStore = zeros(summaryRows, summaryCols, sampleNo);

listing = dir('MEFcutoff*.xlsx');       % prefix for datafiles
fileNumber = size(listing,1);
spreadsheetRowSize = size(S,2)*3;
paramSummary = zeros(fileNumber, spreadsheetRowSize);

% [4] RETRIEVE DATA LOOP

for figCount = 1:size(S,2)  % this loop takes us through cell lines

    column = column+3;          % this is the spreadsheet - specific offset    

    % read in specific data using figureSummary
    % note not all data here are included in the paper
    figureSummaryEXglc = -figureSummary(listing, 'glc_D_t', 'EX_glc_D[e]', column, spreadsheetRowSize, reshapeRow); % '-' to take into account inward EX reactions are -ve
    figureSummaryDMatp = figureSummary(listing, 'atp', 'DM_atp_c_', column, spreadsheetRowSize, reshapeRow);
    figureSummaryDMatp = max(figureSummaryDMatp, 0); % this removes any negative values which might create problems later
    figureSummaryPFK = figureSummary(listing, 'atp', 'PFK', column, spreadsheetRowSize, reshapeRow);
    figureSummaryATPS4mi = figureSummary(listing, 'atp', 'ATPS4mi', column, spreadsheetRowSize, reshapeRow);
    figureSummaryATPyield = figureSummaryDMatp./figureSummaryEXglc;

    % accumulate in 'Store' variables to include all cell lines
    figureSummaryEXglcStore(:,:,figCount)=figureSummaryEXglc; % change sign to avoid confusion of negative lower bound for uptake from exchange reactions
    figureSummaryDMatpStore(:,:,figCount)=figureSummaryDMatp;
    figureSummaryPFKStore(:,:,figCount)=figureSummaryPFK;
    figureSummaryATPS4miStore(:,:,figCount)=figureSummaryATPS4mi;
    figureSummaryATPyieldStore(:,:,figCount)=figureSummaryATPyield;
    
end

% [4] DERIVED PLOTS

% ATP synthesis per molar measured glucose uptake (Qglc values above)
for figCount = 1:sampleNo
    figureSummaryATPyieldMeasuredStore(:,:,figCount)=figureSummaryATPyieldStore(:,:,figCount)*Qglc(figCount);
end

% log ATP synthesis
for figCount = 1:sampleNo
    figureSummarylogATPyieldMeasuredStore(:,:,figCount)=log(figureSummaryATPyieldStore(:,:,figCount)*Qglc(figCount));
end

% glycolysis versus ox phos
for figCount = 1:sampleNo
    figureSummaryGlycolysisStore(:,:,figCount)=figureSummaryPFKStore(:,:,figCount)./figureSummaryATPS4miStore(:,:,figCount);
end


% [5] INITIAL FIGURES

subTitle1 = {};
subTitle1{1} = 'WT';
subTitle1{2} = 'G12D';
subTitle1{3} = 'G12V';
subTitle1{4} = 'Q61L';

mainTitle = 'Estimated maximum ATP capacity';
im = plot4Lines(mainTitle,subTitle1,figureSummaryDMatpStore);
filename = strcat('ATP_capacity4.tif');
saveas(gcf,filename);

mainTitle = 'Estimated exchange of glucose';
im = plot4Lines(mainTitle,subTitle1,figureSummaryEXglcStore);
filename = strcat('Glucose4.tif');
saveas(gcf,filename);

mainTitle = 'Estimated ATP capacity : glucose uptake ratio';
im = plot4Lines(mainTitle,subTitle1,figureSummaryATPyieldStore);
filename = strcat('ATPperGlucose4.tif');
saveas(gcf,filename);

mainTitle = 'Maximum ATP yield from measured glucose uptake (pmol/ng protein/day)';
im = plot4Lines(mainTitle,subTitle1,figureSummaryATPyieldMeasuredStore);
filename = strcat('ATPmeasuredYield4.tif');
saveas(gcf,filename);

%% 
mainTitle = 'log ATP pmol/ng protein/day';
im = plot4Lines(mainTitle,subTitle1,log(figureSummaryATPyieldMeasuredStore));
filename = strcat('lgATPperngPMeasured.tif');
saveas(gcf,filename);

mainTitle = 'PTK / ATPS';
im = plot4Lines(mainTitle,subTitle1,figureSummaryGlycolysisStore);
filename = strcat('glycolysis.tif');
saveas(gcf,filename);


% [6] NEXT FIGURES

subTitle2{1} = 'WT';
subTitle2{2} = 'WT v G12D';
subTitle2{3} = 'WT v G12V';
subTitle2{4} = 'WT v Q61L';

% some derived plots: WT ATP supply and differences

figureSummaryATPdiffsStore(:,:,1) = log(figureSummaryATPyieldMeasuredStore(:,:,1));
figureSummaryATPdiffsStore(:,:,2) = log(figureSummaryATPyieldMeasuredStore(:,:,2))-log(figureSummaryATPyieldMeasuredStore(:,:,1));
figureSummaryATPdiffsStore(:,:,3) = log(figureSummaryATPyieldMeasuredStore(:,:,3))-log(figureSummaryATPyieldMeasuredStore(:,:,1));
figureSummaryATPdiffsStore(:,:,4) = log(figureSummaryATPyieldMeasuredStore(:,:,4))-log(figureSummaryATPyieldMeasuredStore(:,:,1));

mainTitle = 'WT log ATP pmol/ng protein/day and line differences';
im = plot4LinesIndependent(mainTitle,subTitle2,figureSummaryATPdiffsStore);
filename = strcat('lgATPdiffs.tif');
saveas(gcf,filename);


% [7] PROTEIN SYNTHESIS ANALYSIS TAKING INTO ACCOUNT ATP AND GTP

% So take ATPase summary file that includes GTPases as well as ATPases
ATPaseSummary = {};
tRNAlig = zeros(sampleNo,size(Squeeze,2)); % rows = cell lines & columns = different degrees of squeeze
AvoGNo = 6.0221409e+23;     % Avogadro's number
ATP1ngP = 2.19e+13;         % number of ATP's by tRNA ligases for 1ng protein

sheetTitle = {};
sheetTitle{1} = 'WT';
sheetTitle{2} = 'G12D';
sheetTitle{3} = 'G12V';
sheetTitle{4} = 'Q61L';
% extract the fraction of ATPase available for tRNAligases for each squeeze
% value of kcat - corrected expression levels
for sheetCount = 1: sampleNo
    ATPaseSummary = readtable("ATPGTPaseSummary.xlsx", 'Sheet', sheetTitle{sheetCount}); 
    for squeezeCount = 1:size(Squeeze,2)
        % row 24 is tRNA ligases, row 60 is non-ATP/GTP ase abundance
        tRNAlig(sheetCount, squeezeCount) = ATPaseSummary{24,(27+squeezeCount)}/(100-ATPaseSummary{60,(27+squeezeCount)});
    end
end


% now take ATP pmol ng-1.day-1 as basis for calculation (figureSummaryATPyieldMeasuredStore) 

for figCount = 1:sampleNo
    
    sATPyield = figureSummaryATPyieldMeasuredStore(:,:,figCount);   % total ATP available (pmol/ng protein/day)
    stRNAlig = tRNAlig(figCount,:)';                                % fraction for tRNA ligases
    ligATP = repmat(stRNAlig,1,size(sATPyield,2));                  % copy across to align with sATPyield plot
    ligATPc = ligATP.*sATPyield;                                     % map of ATP for ligases (pmol/ng protein/day)
    ligATPm = ligATPc.*(AvoGNo/1e+12);                                % # molecules ATP for ligase /ng protein / day
    ligATPp = ligATPm./ATP1ngP;                                       % ng protein / ng protein /day
    
    figureSummaryPerCentlgATPase(:,:,figCount) = ligATP;            % need to change to fraction not %
    figureSummaryFractionATPlgATPase(:,:,figCount) = ligATPc;
    figureSummaryMolsATPlgATPase(:,:,figCount) = ligATPm;
    figureSummaryProteinRateStore(:,:,figCount) = ligATPp;
end

mainTitle = 'tRNA ligase abundance (#)';
im = plot4Lines(mainTitle,subTitle1,figureSummaryPerCentlgATPase);
filename = strcat('tRNAligaseAbundance.tif');
saveas(gcf,filename);

mainTitle = 'ATP available for tRNA ligase';
im = plot4Lines(mainTitle,subTitle1,figureSummaryFractionATPlgATPase);
filename = strcat('ATPfortRNAligase.tif');
saveas(gcf,filename);

mainTitle = 'Protein/day';
im = plot4Lines(mainTitle,subTitle1,figureSummaryProteinRateStore);
filename = strcat('ProteinPerDay.tif');
saveas(gcf,filename);


% Further derived plot: WT Protein and differences

figureSummaryProteindiffsStore(:,:,1) = figureSummaryProteinRateStore(:,:,1);
figureSummaryProteindiffsStore(:,:,2) = figureSummaryProteinRateStore(:,:,2)-figureSummaryProteinRateStore(:,:,1);
figureSummaryProteindiffsStore(:,:,3) = figureSummaryProteinRateStore(:,:,3)-figureSummaryProteinRateStore(:,:,1);
figureSummaryProteindiffsStore(:,:,4) = figureSummaryProteinRateStore(:,:,4)-figureSummaryProteinRateStore(:,:,1);

mainTitle = 'WT Protein production rate and differences for mutant lines';
im = plot4LinesIndependent(mainTitle,subTitle2,figureSummaryProteindiffsStore);
filename = strcat('Proteindiffs.tif');
saveas(gcf,filename);

% From Figure 1G modified as mean 24-48 and 48-72 hours ug/ml per ug/ml per
% day (ie a rate)

% Qpr = [0.722	1.39	1.71	0.889]'; Not sure these are correct
Qpr = [0.326	.979	.651	0.741]';

r2 = zeros(size(figureSummaryProteinRateStore, 1), size(figureSummaryProteinRateStore, 2));
for rowNo = 1:size(r2,1)
    for colNo = 1:size(r2,2)
        fit = fitlm(Qpr, squeeze(figureSummaryProteinRateStore(rowNo,colNo,:)));
        r2(rowNo, colNo)=fit.Rsquared.Ordinary;
    end
end

clf
image(r2,'CDataMapping','scaled');
title('r2 for relationship between measured and predicted protein rates');
xlabel('EX_ v protein abundance');
ylabel('Compression');
xticks([ 4  6  8  10 ]);
xticklabels([0.025  0.1  0.4  1.6 ]);
yticks([2 4 6 8 10 12]);
yticklabels([0.0005 0.002 0.008 0.032 0.128 0.5]);
colorbar
filename = strcat('ProteinRatesR2.tif');
saveas(gcf,filename);

meanProteinRate = squeeze(mean(figureSummaryProteinRateStore,3));
meanQpr = mean(Qpr);
turnover = 0.5; % fraction of cellular proteins turned over per day
totalQpr = meanQpr + turnover; % https://bionumbers.hms.harvard.edu/bionumber.aspx?id=109937

for figCount = 1:sampleNo
    climMax = max(figureSummaryProteinRateStore, [], 'all');
    clim = [0, climMax];   
    subplot(2,2,figCount);
    im = image(figureSummaryProteinRateStore(:,:,figCount),'CDataMapping','scaled');
    sgtitle('Estimated protein rates with measured rate contour');   
    % title(subTitle(count));
    xlabel('EX_ v protein abundance');
    ylabel('Compression');
    xticks([ 4  6  8  10 ]);
    xticklabels([0.025  0.1  0.4  1.6 ]);
    yticks([2 4 6 8 10 12]);
    yticklabels([0.0005 0.002 0.008 0.032 0.128 0.5]);
    colorbar
    % caxis(clim);
    pbaspect([1 1 1]);   
    hold on
    contour(figureSummaryProteinRateStore(:,:,figCount), [Qpr(figCount) Qpr(figCount)], 'LineWidth', 2, 'Color', 'red' );
    hold off   
end    

filename = strcat('ProteinRateContour.tif');
saveas(gcf,filename);

clf
im=image(meanProteinRate,'CDataMapping','scaled');
title('Mean Protein Turnover Rate');
xlabel('EX_ v protein abundance');
ylabel('Compression');
xticks([ 4  6  8  10 ]);
xticklabels([0.025  0.1  0.4  1.6 ]);
yticks([2 4 6 8 10 12]);
yticklabels([0.0005 0.002 0.008 0.032 0.128 0.5]);
hold on
contour(meanProteinRate, [meanQpr meanQpr], 'LineWidth', 2, 'Color', 'red' );
colorbar
hold off

filename = strcat('meanProteinRateContour.tif');
saveas(gcf,filename);
