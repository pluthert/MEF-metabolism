function figureSummaryData = figureSummary(listing, sheetName, rowName, column, spreadsheetRowSize, reshapeRow)

%   fileNumber = number of files in directory to be analysed
%   listing = directory listing of relevant files in directory
%   sheetName = spreadsheet name
%   rowName = row of interest in sheet sheetName
%   column = column of spreadsheet to be read
%   spreadsheetRowSzie = width of spreadsheet row

fileNumber = size(listing,1);
paramSummary = zeros(fileNumber,spreadsheetRowSize);

    for i=1:fileNumber
        file = listing(i).name;
        tab = readtable(file, 'sheet', sheetName, 'ReadRowNames', true, 'ReadVariableNames', true);
            if ~isempty(tab{:,1:end-1})       % deals with situation where no solution found (tab empty)
                if max(strcmp(rowName,tab.Properties.RowNames))
                    paramRow = tab{rowName,1:end-1};
                else
                    paramRow = zeros(1, spreadsheetRowSize); 
                end    
            else
                paramRow = zeros(1, spreadsheetRowSize); 
            end
        paramSummary(i,:) = paramRow;
    end
    figureSummaryData = reshape(paramSummary(:,column), [reshapeRow],[]); 

end

