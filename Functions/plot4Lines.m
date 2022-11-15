function im = plot4Lines(mainTitle,subTitle,dataSummary)
% mainTitle = string for overarching title
% subTitle = string array (4) for each plot
% dataSummary = figureSummaryStore dim1 = rows, dim2 = cols, dim3 = lines
% for now graphing parameters 'hard' coded

for count = 1:4
    climMax = max(dataSummary, [], 'all');
    clim = [0, climMax];
    subplot(2,2,count);
    im=image(dataSummary(:,:,count),'CDataMapping','scaled');
    sgtitle(mainTitle);
    title(subTitle(count));
    xlabel('EX_ v protein abundance');
    ylabel('Compression');
    xticks([ 4  6  8  10 ]);
    xticklabels([0.025  0.1  0.4  1.6 ]);
    yticks([2 4 6 8 10 12]);
    yticklabels([0.0005 0.002 0.008 0.032 0.128 0.5]);
    colorbar
    caxis(clim);
    pbaspect([1 1 1]);
    
end

end

