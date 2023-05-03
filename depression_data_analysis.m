clear
close all
clc

%% load score
basePath = pwd;
addpath(genpath(fullfile(basePath, 'Dataset', 'DepressionDataset')))
score = readtable("scores.csv");
[nSubj, nFeat] = size(score);

% conditions: 1 - 23
% controls: 24 - 55

%% replace NaN in madrs of controls (madrs < 10 = absence of depression symptoms)
rng(42)
score(24:end, end - 1) = cellstr(string(randi([0 9],32,1)));
score(24:end, end) = cellstr(string(randi([0 9],32,1)));

% SBAGLIATO

%% create the matrix X containing: days, gender, age, afftype, melanch, inpatient, edu, marrage, work, madrs1, madrs2
X = nan(nSubj, nFeat - 1); % we don't consider the covariate "number"

for i = [2, 3]
    X(:, i - 1) = table2array(score(:,i));
end

for i = [4, 8]
    
    slot = string(unique(score{:,i}));
    slot(slot == "") = [];

    for j = 1:length(slot)
        X(string(score{:,i}) == slot(j), i - 1) = j;
    end 

end

for i = [5, 6, 7, 9, 10, 11, 12]
    X(:, i - 1) = str2double(score{:, i});
end

% subsitute madrs1 and madrs2 with their mean
madrs = round(mean([X(:,end -1) X(:,end)]'))';
X(:,[end - 1, end]) = [];
X = [X madrs];

clear i j madrs nFeat nSubj slot
%% distibution of gender and age
figure
% gender
subplot(121)
b = bar([sum(X(1:23,2) == 1) sum(X(24:end,2) == 1); sum(X(1:23,2) == 2) sum(X(24:end,2) == 2)]);
ylim([0 25])
xticks(1:2)
xticklabels({'female', 'male'})
legend('depression', 'normal', 'Location', 'northwest')
title('gender')
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% age
subplot(122)
x = [];
for i = 1:10
    x = [x; sum(X(1:23,3) == i) sum(X(24:end,3) == i)];
end
b = bar(x);
ylim([0 7])
xticks(1:10)
xticklabels(unique(score{:,4}))
legend('depression', 'normal', 'Location', 'northwest')
title('age')
for i = 1:length(b)
    xtips = b(i).XEndPoints;
    ytips = b(i).YEndPoints;
    labels = string(b(i).YData);
    text(xtips,ytips,labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end

clear b i labels labels1 labels2 x xtips xtips1 xtips2 ytips2 ytips1 ytips
%% distibution of afftype, melanch, inpatient, edu, marriage, work, madrs
figure
% afftype
subplot(241)
b = bar([sum(X(1:23,4) == 1) sum(X(1:23,4) == 2) sum(X(1:23,4) == 3)]);
ylim([0 25])
xticks(1:3)
xticklabels({'bipolar II', 'unipolar depressive', 'bipolar I'})
xtickangle(45)
legend('depression', 'Location', 'northwest')
title('type of disorder')
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips, ytips, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% melanch
subplot(242)
b = bar([sum(X(1:23,5) == 1) sum(X(1:23,5) == 2) sum(isnan(X(1:23,5)))]);
ylim([0 25])
xticks(1:3)
xticklabels({'yes', 'no', 'NaN'})
xtickangle(45)
legend('depression', 'Location', 'northwest')
title('presence of melancholia')
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips, ytips, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% inpatient
subplot(243)
b = bar([sum(X(1:23,6) == 1) sum(X(1:23,6) == 2)]);
ylim([0 25])
xticks(1:2)
xticklabels({'yes', 'no'})
xtickangle(45)
legend('depression', 'Location', 'northwest')
title('hospitalization')
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips, ytips, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% edu
subplot(244)
b = bar([sum(X(1:23,7) == 1) sum(X(1:23,7) == 2) sum(X(1:23,7) == 3) sum(isnan(X(1:23,7)))]);
ylim([0 25])
xticks(1:4)
x = unique(score{1:23,8});
x(end + 1) = cellstr('NaN');
xticklabels(x(2:end))
xtickangle(45)
legend('depression', 'Location', 'northwest')
title('education')
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips, ytips, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% marriage
subplot(245)
b = bar([sum(X(1:23,7) == 1) sum(X(1:23,7) == 2)]);
ylim([0 25])
xticks(1:2)
xticklabels({'married', 'single'})
xtickangle(45)
legend('depression', 'Location', 'northwest')
title('marriage')
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips, ytips, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% work
subplot(246)
b = bar([sum(X(1:23,8) == 1) sum(X(1:23,8) == 2)]);
ylim([0 25])
xticks(1:2)
xticklabels({'working', 'unemployed'})
xtickangle(45)
legend('depression', 'Location', 'northwest')
title('work')
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips, ytips, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% madrs
countMadrs = [];
slot = unique(X(1:23, end));
for i = 1:length(slot)
    countMadrs = [countMadrs sum(X(1:23, end) == slot(i))];
end
subplot(247)
b = bar(countMadrs);
ylim([0 25])
xticks(1:length(slot))
xticklabels(cellstr(string(slot)))
xtickangle(45)
legend('depression', 'Location', 'northwest')
title('MADRS score') % this MADRS score has been obtained computing the average
                     % between the MADRS score when the measurement started and
                     % the MADRS score when the measurement stopped
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(b.YData);
text(xtips, ytips, labels,'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

clear b countMadrs i labels slot x xtips ytips

%% load conditions
cndTables = dir(fullfile(basePath, 'Dataset', 'DepressionDataset', 'condition'));
cndTables = cndTables(3:end);

sp = split(string({cndTables.name}.'), ["_","."]);
num = str2double(sp(:,2));
[~, idx] = sort(num);
cndTables = cndTables(idx);

clear sp num idx
%% heatmap conditions

% we want to visualize the motor activity across days of week and hours of
% day

T = [];

for cnd = 1:length(cndTables) % for each subjects
    
    % load data
    currCnd = readtable(cndTables(cnd).name);

    disp(['working on ', cndTables(cnd).name])

    % add to currCnd the columns hour, dayNumber, dayName
    timeSamp = string(currCnd{:,1}); timeSamp = split(timeSamp, " ");
    data = timeSamp(:,1);
    [dayNumber, dayName] = weekday(data);
    hour = timeSamp(:,2); hour = split(hour, ":"); hour = str2double(hour(:,1));
    currCnd = addvars(currCnd, dayNumber, dayName, hour);
    
     T = [T; currCnd(:,3:6)];

end     

T = sortrows(T,[2 4]);
figure
heatmap(T, "dayNumber", "hour", "ColorMethod", "mean", "ColorVariable","activity")

clear cnd currCnd data dayName dayNumber hour timeSamp T
%% load controls
ctrTables = dir(fullfile(basePath, 'Dataset', 'DepressionDataset', 'control'));
ctrTables = ctrTables(3:end);

sp = split(string({ctrTables.name}.'), ["_","."]);
num = str2double(sp(:,2));
[~, idx] = sort(num);
ctrTables = ctrTables(idx);

clear sp num idx
%% heatmap controls
T = [];

for ctr = 1:length(ctrTables) % for each subjects
    
    % load data
    currCtr = readtable(ctrTables(ctr).name);

    disp(['working on ', ctrTables(ctr).name])

    % add to currCtr the columns hour, dayNumber, dayName
    timeSamp = string(currCtr{:,1}); timeSamp = split(timeSamp, " ");
    data = timeSamp(:,1);
    [dayNumber, dayName] = weekday(data);
    hour = timeSamp(:,2); hour = split(hour, ":"); hour = str2double(hour(:,1));
    currCtr = addvars(currCtr, dayNumber, dayName, hour);

    T = [T; currCtr(:,3:6)];

end 

T = sortrows(T,[2 4]);
figure
heatmap(T, "dayNumber", "hour", "ColorMethod", "mean", "ColorVariable","activity")

clear ctr currCtr data dayName dayNumber hour timeSamp T
%% 1a) Is there any difference between controls and depressed patients in daily activity (or any derived score)?

%% decide which days we want to keep

% condition
actualDaysCnd = [];
for cnd = 1:length(cndTables) % for each subjects
    
    % load data
    currCnd = readtable(cndTables(cnd).name);

    disp(['working on ', cndTables(cnd).name])
    disp(['days: ', num2str(score{cnd,2})])

    % add to currCnd the columns hour, dayNumber, dayName
    timeSamp = string(currCnd{:,1}); timeSamp = split(timeSamp, " ");
    data = timeSamp(:,1);
    [dayNumber, dayName] = weekday(data);
    hour = timeSamp(:,2); hour = split(hour, ":"); hour = str2double(hour(:,1));
    currCnd = addvars(currCnd, dayNumber, dayName, hour);

    days = 1;
    i = 1; % counts for rows of currCnd

    while i < size(currCnd, 1) % for each day of recording
        
        idx_i = find(currCnd{:,2} == currCnd{i,2}); 
        day = currCnd(idx_i,:);

        ii = 1; % counts for rows of day
        sumHour = [];
        hours = 1; % count for the hours (index for sumHour)
       
        while ii < size(day,1) % for each hour in a day of recording

            timeSamp = string(day{ii,1}); timeSamp = split(timeSamp, " ");
            h = split(timeSamp(2), ":"); h = str2double(h(1));
    
            idx_ii = find(day{:,6} == h);
            sumHour(hours) = sum(day{idx_ii, 3});
            
            ii = ii + length(idx_ii); % go to the next hour
            hours = hours + 1;

        end % while ii

        if sum(sumHour > 0) <= 10 % if the device has been wore for less or equal than 10 hours
            % the samples of the current day must be discarded
            disp(['day to be discarded: ', num2str(days)])
        end % if

        i = i + length(idx_i); % go to the next day
        days = days + 1;

    end % while i

    actualDaysCnd(cnd) = days - 1;
    disp(['actual days: ', num2str(actualDaysCnd(cnd))])
    disp(" ")

end % for cnd

% controls
actualDaysCtr = [];
for ctr = 1:length(ctrTables) % for each subjects
    
    % load data
    currCtr = readtable(ctrTables(ctr).name);

    disp(['working on ', ctrTables(ctr).name])
    disp(['days: ', num2str(score{19 + ctr,2})])

    % add to currCtr the columns hour, dayNumber, dayName
    timeSamp = string(currCtr{:,1}); timeSamp = split(timeSamp, " ");
    data = timeSamp(:,1);
    [dayNumber, dayName] = weekday(data);
    hour = timeSamp(:,2); hour = split(hour, ":"); hour = str2double(hour(:,1));
    currCtr = addvars(currCtr, dayNumber, dayName, hour);

    days = 1;
    i = 1; % counts for rows of currCtr

    while i < size(currCtr, 1) % for each day of recording
        
        idx_i = find(currCtr{:,2} == currCtr{i,2}); 
        day = currCtr(idx_i,:);

        ii = 1; % counts for rows of day
        sumHour = [];
        hours = 1; % count for the hours (index for sumHour)
       
        while ii < size(day,1) % for each hour in a day of recording

            timeSamp = string(day{ii,1}); timeSamp = split(timeSamp, " ");
            h = split(timeSamp(2), ":"); h = str2double(h(1));
    
            idx_ii = find(day{:,6} == h);
            sumHour(hours) = sum(day{idx_ii, 3});
            
            ii = ii + length(idx_ii); % go to the next hour
            hours = hours + 1;

        end % while ii

        if sum(sumHour > 0) <= 10 % if the device has been worn for less or equal than 10 hours
            % the samples of the current day must be discarded
            disp(['day to be discarded: ', num2str(days)])
        end % if

        i = i + length(idx_i); % go to the next day
        days = days + 1;

    end % while i

    actualDaysCtr(ctr) = days - 1;
    disp(['actual days: ', num2str(actualDaysCtr(ctr))])
    disp(" ")

end % for ctr

% looking at which were the days to be discarded according to our criterion
% we saw that these days were for the most of the subjects the first and
% the lasts
% we decide to discard the first days if marked as a day to be discarded, 
% all the days after the first of the lasts days marked as days to be
% discarded and the the first of the lasts days marked as days to be
% discarded because doing this the remaining days were enough for all the
% subjects!

clear actualDaysCnd actualDaysCtr cnd ctr currCnd currCtr data day dayName dayNumber days h hours hour i idx_i idx_ii ...
    ii sumHour timeSamp
%% calculate biomarkers for CONDITIONS discarding days
meanCnd = nan(1, length(cndTables));
sdCnd = nan(1, length(cndTables));
cvCnd = nan(1, length(cndTables));
weekdayMeanCnd = nan(1, length(cndTables));
weekdaySdCnd = nan(1, length(cndTables));
weekdayCvCnd = nan(1, length(cndTables));
weekendMeanCnd = nan(1, length(cndTables));
weekendSdCnd = nan(1, length(cndTables));
weekendCvCnd = nan(1, length(cndTables));
autocorrDailyMeanCnd = zeros(1, length(cndTables));

for cnd = 1:length(cndTables) % for each subjects
    
    % load data
    currCnd = readtable(cndTables(cnd).name);

    disp(['working on ', cndTables(cnd).name])

    % add to currCnd the columns hour, dayNumber, dayName
    timeSamp = string(currCnd{:,1}); timeSamp = split(timeSamp, " ");
    data = timeSamp(:,1);
    [dayNumber, dayName] = weekday(data);
    hour = timeSamp(:,2); hour = split(hour, ":"); hour = str2double(hour(:,1));
    currCnd = addvars(currCnd, dayNumber, dayName, hour);

    %% discarding days
    days = 1;
    i = 1; % counts for rows of currCnd
    idxFirstDay = [];
    countFirstDay = 0;
    while i < size(currCnd, 1) % for each day of recording
        
        idx_i = find(currCnd{:,2} == currCnd{i,2}); 
        day = currCnd(idx_i,:);

        ii = 1; % counts for rows of day
        sumHour = [];
        hours = 1; % count for the hours (index for sumHour)
       
        while ii < size(day,1) % for each hour in a day of recording

            timeSamp = string(day{ii,1}); timeSamp = split(timeSamp, " ");
            h = split(timeSamp(2), ":"); h = str2double(h(1));
    
            idx_ii = find(day{:,6} == h);
            sumHour(hours) = sum(day{idx_ii, 3});
            
            ii = ii + length(idx_ii); % go to the next hour
            hours = hours + 1;

        end % while ii

        if sum(sumHour > 0) <= 10 &&  days == 1 % if the device has been worn for less or equal than 10 hours and the current day is the first one
            % the samples of the current day must be discarded
            idxFirstDay = idx_i;
            countFirstDay = 1;
        elseif sum(sumHour > 0) <= 10 &&  days ~= 1  % if the device has been worn for less or equal than 10 hours and the current day is not the first one
            % the samples of the current day and of all the following days must be discarded
            firstIdxToDelete = idx_i(1);
            break
        end % if

        i = i + length(idx_i); % go to the next day
        days = days + 1;

    end % while i

    currCnd([idxFirstDay; (firstIdxToDelete:size(currCnd, 1))'],:) = [];
    actualDaysCnd(cnd) = days - 1 - countFirstDay;

    %% calculate biomarkers
    i = 1; % counts for rows of currCnd
    dailyMean = [];
    days = 1; % index of dailyMean
    
    while i < size(currCnd, 1) % for each day of recording

        idx_i = find(currCnd{:,2} == currCnd{i,2}); 
        day = currCnd(idx_i,:);
        dailyMean(days) = mean(day{:,3});
        
        i = i + length(idx_i); % go to the next day
        days = days + 1;
    
    end % while i

    [acf, lags, bounds] = autocorr(dailyMean);
    if any(abs(acf(2:end)) > bounds(1))
        autocorrDailyMeanCnd(cnd) = 1; % the past explains for the future
    end 

    meanCnd(cnd) = mean(currCnd{:,3});
    sdCnd(cnd) = std(currCnd{:,3});
    cvCnd(cnd) = sdCnd(cnd)/abs(meanCnd(cnd));

    idxWeekend = currCnd{:,"dayNumber"} == 7 | currCnd{:,"dayNumber"} == 1; % saturday = 7; sunday = 1
    
    weekdayMeanCnd(cnd) = mean(currCnd{~idxWeekend,3});
    weekdaySdCnd(cnd) = std(currCnd{~idxWeekend,3});
    weekdayCvCnd(cnd) = weekdaySdCnd(cnd)/abs(weekdayMeanCnd(cnd));

    weekendMeanCnd(cnd) = mean(currCnd{idxWeekend,3});
    weekendSdCnd(cnd) = std(currCnd{idxWeekend,3});
    weekendCvCnd(cnd) = weekendSdCnd(cnd)/abs(weekendMeanCnd(cnd));

end % for cnd

%clear acf bounds cnd countFirstDay currCnd dailyMean data day dayName dayNumber days firstIdxToDelete ...
%    h hours hour i idx_ii idx_i idxFirstDay idxWeekend ii lags sumHour timeSamp 
%% calculate biomarkers for CONTROLS discarding days
meanCtr = nan(1, length(ctrTables));
sdCtr = nan(1, length(ctrTables));
cvCtr = nan(1, length(ctrTables));
weekdayMeanCtr = nan(1, length(ctrTables));
weekdaySdCtr = nan(1, length(ctrTables));
weekdayCvCtr = nan(1, length(ctrTables));
weekendMeanCtr = nan(1, length(ctrTables));
weekendSdCtr = nan(1, length(ctrTables));
weekendCvCtr = nan(1, length(ctrTables));
autocorrDailyMeanCtr = zeros(1, length(ctrTables));

for ctr = 1:length(ctrTables) % for each subjects
    
    % load data
    currCtr = readtable(ctrTables(ctr).name);

    disp(['working on ', ctrTables(ctr).name])

    % add to currCtr the columns hour, dayNumber, dayName
    timeSamp = string(currCtr{:,1}); timeSamp = split(timeSamp, " ");
    data = timeSamp(:,1);
    [dayNumber, dayName] = weekday(data);
    hour = timeSamp(:,2); hour = split(hour, ":"); hour = str2double(hour(:,1));
    currCtr = addvars(currCtr, dayNumber, dayName, hour);

    %% discarding days
    days = 1;
    i = 1; % counts for rows of currCtr
    idxFirstDay = [];
    countFirstDay = 0;
    while i < size(currCtr, 1) % for each day of recording
        
        idx_i = find(currCtr{:,2} == currCtr{i,2}); 
        day = currCtr(idx_i,:);

        ii = 1; % counts for rows of day
        sumHour = [];
        hours = 1; % count for the hours (index for sumHour)
       
        while ii < size(day,1) % for each hour in a day of recording

            timeSamp = string(day{ii,1}); timeSamp = split(timeSamp, " ");
            h = split(timeSamp(2), ":"); h = str2double(h(1));
    
            idx_ii = find(day{:,6} == h);
            sumHour(hours) = sum(day{idx_ii, 3});
            
            ii = ii + length(idx_ii); % go to the next hour
            hours = hours + 1;

        end % while ii

        if sum(sumHour > 0) <= 10 &&  days == 1 % if the device has been worn for less or equal than 10 hours and the current day is the first one
            % the samples of the current day must be discarded
            idxFirstDay = idx_i;
            countFirstDay = 1;
        elseif sum(sumHour > 0) <= 10 &&  days ~= 1  % if the device has been worn for less or equal than 10 hours and the current day is not the first one
            % the samples of the current day and of all the following days must be discarded
            firstIdxToDelete = idx_i(1);
            break
        end % if

        i = i + length(idx_i); % go to the next day
        days = days + 1;

    end % while i

    currCtr([idxFirstDay; (firstIdxToDelete:size(currCtr, 1))'],:) = [];
    actualDaysCtr(ctr) = days - 1 - countFirstDay;

    %% computing biomarkers
    i = 1; % counts for rows of currCtr
    dailyMean = [];
    days = 1; % index of dailyMean
    
    while i < size(currCtr, 1) % for each day of recording

        idx_i = find(currCtr{:,2} == currCtr{i,2}); 
        day = currCtr(idx_i,:);
        dailyMean(days) = mean(day{:,3});
        
        i = i + length(idx_i); % go to the next day
        days = days + 1;
    
    end % while i

    [acf, lags, bounds] = autocorr(dailyMean);
    if any(abs(acf(2:end)) > bounds(1))
        autocorrDailyMeanCtr(ctr) = 1; % the past explains for the future
    end
    
    meanCtr(ctr) = mean(currCtr{:,3});
    sdCtr(ctr) = std(currCtr{:,3});
    cvCtr(ctr) = sdCtr(ctr)/abs(meanCtr(ctr));

    idxWeekend = currCtr{:,"dayNumber"} == 7 | currCtr{:,"dayNumber"} == 1; % saturday = 7; sunday = 1
    
    weekdayMeanCtr(ctr) = mean(currCtr{~idxWeekend,3});
    weekdaySdCtr(ctr) = std(currCtr{~idxWeekend,3});
    weekdayCvCtr(ctr) = weekdaySdCtr(ctr)/abs(weekdayMeanCtr(ctr));

    weekendMeanCtr(ctr) = mean(currCtr{idxWeekend,3});
    weekendSdCtr(ctr) = std(currCtr{idxWeekend,3});
    weekendCvCtr(ctr) = weekendSdCtr(ctr)/abs(weekendMeanCtr(ctr));


end % for ctr

clear acf bounds ctr countFirstDay currCtr dailyMean data day dayName dayNumber days firstIdxToDelete ...
    h hours hour i idx_ii idx_i idxFirstDay idxWeekend ii lags sumHour timeSamp 
%% uploade the column days in the table score and in the matrix X

score(:,2) = array2table([actualDaysCnd actualDaysCtr]');
X(:,1) = [actualDaysCnd actualDaysCtr]';

clear actualDaysCtr actualDaysCnd
%%

Mean = [meanCnd, meanCtr];
sd = [sdCnd, sdCtr];
cv= [cvCnd, cvCtr];
weekdayMean = [weekdayMeanCnd, weekdayMeanCtr];
weekdaySd = [weekdaySdCnd, weekdaySdCtr];
weekdayCv = [weekdayCvCnd, weekdayCvCtr];
weekendMean = [weekendMeanCnd, weekendMeanCtr];
weekendSd = [weekendSdCnd, weekendSdCtr];
weekendCv = [weekendCvCnd, weekendCvCtr];
autocorrDailyMean =  [autocorrDailyMeanCnd, autocorrDailyMeanCtr]; 

clear meanCtr meanCnd sdCtr sdCnd cvCnd cvCtr weekdayMeanCtr weekdayMeanCnd weekdaySdCnd ...
    weekdaySdCtr weekdayCvCnd weekdayCvCtr weekendMeanCnd weekendMeanCtr weekendSdCtr weekendSdCnd ...
    weekendCvCtr weekendCvCnd autocorrDailyMeanCnd autocorrDailyMeanCtr
%% boxplot

% BOXPLOT works with grouping variables, so you can manually append all
% of your data together and then create a grouping variable that lets
% boxplot know which belongs to first and which for second
potentialBiomarkers = [Mean; sd; cv; weekdayMean; weekdaySd; weekdayCv; weekendMean; weekendSd; weekendCv];
potentialBiomarkersNames = ["mean", "sd", "cv", "weekdayMean", "weekdaySd", "weekdayCv", "weekendMean", "weekendSd", "weekendCv"];
grp = [zeros(1,length(cndTables)),ones(1,length(ctrTables))];
figure
for i = 1:size(potentialBiomarkers,1)
    subplot(3,3,i)
    boxplot(potentialBiomarkers(i,:), grp, 'Labels', {'depression','normal'})
    title(potentialBiomarkersNames(i))
end
sgtitle(' biomarkers differences across states (depression/normal)')

% by visual inspection we can see differences, but we have to perform
% statistical tests to verify if this differences are statisticative

figure
b = bar([sum(autocorrDailyMean(1:23) == 0) sum(autocorrDailyMean(1:23) == 1); sum(autocorrDailyMean(24:end) == 0) sum(autocorrDailyMean(24:end) == 1)]);
ylim([0 27])
xticks(1:2)
xticklabels({'depression', 'normal'})
legend('no autocorr', 'autocorr', 'Location', 'northwest')
title('autocorr(mean daily activity)')
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

trueLabel = [zeros(1,length(cndTables)),ones(1,length(ctrTables))];
predLabel = autocorrDailyMean;
cm = confusionmat(trueLabel, predLabel);
accuracy = sum(diag(cm))/sum(cm(:));

% from this plot and from the computed accuracy we can see that autocorrelation would not be helpful for
% distinguishing between conditions and controls, thus we will not
% consider this variable anymore

clear b grp i labels1 labels2 xtips1 xtips2 ytips1 ytips2
%% statistical tests

h = nan(1, size(potentialBiomarkers,1));
p = nan(1, size(potentialBiomarkers,1));

for i = 1:size(potentialBiomarkers,1)
    if ~ lillietest(potentialBiomarkers(i,1:23)) && ~lillietest(potentialBiomarkers(i,24:end)) % if the variables are normally distributed
        [h(i), p(i)] = ttest2(potentialBiomarkers(i,1:23), potentialBiomarkers(i,24:end));
    else
        [p(i), h(i)] = ranksum(potentialBiomarkers(i,1:23), potentialBiomarkers(i,24:end));
    end
end

[h, ~, ~, p] = fdr_bh(p, 0.05,'pdep','on');

% there are significant differences betweek depression and normal state for
% each biomarkers

clear i
%% boxplot + p-values

grp = [zeros(1,length(cndTables)),ones(1,length(ctrTables))];
figure
for i = 1:size(potentialBiomarkers,1)
    subplot(3,3,i)
    boxplot(potentialBiomarkers(i,:), grp, 'Labels', {'depression','normal'})
    title(potentialBiomarkersNames(i))

    if i == 1 || i == 4 || i == 7
        ylim([40 450])
        text(1.3, 400, ['P = ', num2str(round(p(i),3))])
    elseif i == 2 || i == 5 || i == 8
        ylim([110 680])
        text(1.3, 600, ['P = ', num2str(round(p(i),3))])
    elseif i == 3 || i == 6 || i == 9
        ylim([1.1 3])
        text(1.3, 2.8, ['P = ', num2str(round(p(i),3))])
    end

end
sgtitle(' biomarkers differences across states (depression/normal)')

clear grp i
%% 1b) Are the motor activity or any derived scores associated to depression score severity?

r = nan(1, size(potentialBiomarkers,1));
pCorr = nan(1, size(potentialBiomarkers,1));

for i = 1:size(potentialBiomarkers,1)
    if ~lillietest(potentialBiomarkers(i,1:23)) && ~lillietest(potentialBiomarkers(i,24:end)) % if the variables are normally distributed
        [r(i), pCorr(i)] = corr(potentialBiomarkers(i,:)', X(:,end)); % corr compute the Pearson's correlations between columns
    else
        [r(i), pCorr(i)] = corr(potentialBiomarkers(i,:)', X(:,end), 'Type', 'Spearman'); % corr compute the Spearman's correlations between columns
    end
end

[hCorr, ~, ~, pCorr] = fdr_bh(pCorr, 0.05,'pdep','on');

% r measures type and magnitude of the linear relation between two
% variable

% p < 0.5 => h = 1 => correlation

% Mean, sd, weekdayMean, weekdaySd, weekendMean, weekendSd are correlated with madrs

%% plot

figure
j = 1;
for i = [1 2 4 5 7 8]
    subplot(3,2,j)
    plot(X(:,end), potentialBiomarkers(i,:), 'o', 'LineWidth', 2)
    grid on
    hold on
    mod = fitlm(X(:,end), potentialBiomarkers(i,:));
    y = mod.Coefficients{1,1} + mod.Coefficients{2,1}*X(:,end);
    plot(X(:,end), y, 'r', 'LineWidth', 2)

    ylim([0 700])
    text(10, 620, ['rho = ',num2str(round(r(i),3))])
    text(10, 550, ['p = ', num2str(round(pCorr(i),3))])

    xlabel('madrs')
    ylabel(potentialBiomarkersNames(i))
    hold off
    j = j + 1;
    clear mod y
end

sgtitle('significant correlations between biomarkers and madrs')

clear i j biom mod
%% 2) Can the daily activity (or any derived score) capable of distinguishing patients from controls?

auc = nan(1,size(potentialBiomarkers,1));
legend_label = cell(1,size(potentialBiomarkers,1) + 1);

figure
for i = 1:size(potentialBiomarkers,1)

    % frequency distribution
    bin_size = .25;
    [x0_axis, x0_freq] = compute_prob_density(potentialBiomarkers(i,1:19), bin_size); % or replace it with histogram function
    [x1_axis, x1_freq] = compute_prob_density(potentialBiomarkers(i,20:end), bin_size);
    
    % compute ROC
    [fpr, tpr] = compute_roc_curve(x0_axis, x0_freq, x1_axis, x1_freq, bin_size);
    auc(i) = trapz(fpr, tpr);
    if auc(i) < .5
        auc(i) = 1 - auc(i);
        tmp = fpr;
        fpr = tpr;
        tpr = tmp;
    end
   
    plot(fpr, tpr, '.-')
    hold on
    legend_label{i} = [potentialBiomarkersNames{i}, ', auc: ', num2str(auc(i))];
   
end

plot([0,1], [0,1], '--k')
legend_label{i+1} = 'chance level';
legend(legend_label, 'Location', 'southeast')
xlabel('fpr (1-specificity)')
ylabel('tpr (sensitivity)')
xlim([-.1, 1.1])
ylim([-.1, 1.1])

% based on auc metric the biomarkers seems to distinguish between conditions and controls
% and the weekend metric seems to be the better

clear bin_size fpr tpr i legend_label tmp x0_axis x0_freq x1_axis x1_freq
%% sensitivity analysis

rCovar = nan(3,size(potentialBiomarkers,1));
pCovar = nan(3,size(potentialBiomarkers,1));
hCovar = nan(3,size(potentialBiomarkers,1));

for i = 1:3 % days, gender, age
   
    for j = 1:size(potentialBiomarkers,1)
        if ~lillietest(potentialBiomarkers(j,1:23)) && ~lillietest(potentialBiomarkers(j,24:end)) % if the variables are normally distributed
            [rCovar(i,j), pCovar(i,j)] = corr(potentialBiomarkers(j,:)', X(:,i)); % corr compute the Pearson's correlations between columns
        else
            [rCovar(i,j), pCovar(i,j)] = corr(potentialBiomarkers(j,:)', X(:,i), 'Type', 'Spearman'); % corr compute the Spearman's correlations between columns
        end
        
    end
    % here we are testing the correlation of the current covariate against all biomarkers
    [hCovar(i,:), ~, ~, pCovar(i,:)] = fdr_bh(pCovar(i,:), 0.05, 'pdep', 'on');
end

% sd and weekdaysSd are significantly correlated with the variable age

% we have considered only the variables without missing data for both the depression
% and the normal state because we want to draw conclusion considering both
% the states

clear i j
%% plot

figure
% sd vs. age
subplot(121)
plot(X(:,3), potentialBiomarkers(2,:), 'o', 'LineWidth', 2)
grid on
hold on
mod = fitlm(X(:,3), potentialBiomarkers(2,:));
y = mod.Coefficients{1,1} + mod.Coefficients{2,1}*X(:,3);
plot(X(:,3), y, 'r', 'LineWidth', 2)
ylim([0 700])
text(4, 620, ['rho = ',num2str(round(rCovar(3,2),3))])
text(4, 580, ['p = ', num2str(round(pCovar(3,2),3))])
xlabel('age')
xticks(1:10)
xticklabels(unique(score{:,4}))
ylabel('sd')

% weekdaySd vs. age
subplot(122)
plot(X(:,3), potentialBiomarkers(5,:), 'o', 'LineWidth', 2)
grid on
hold on
mod = fitlm(X(:,3), potentialBiomarkers(5,:));
y = mod.Coefficients{1,1} + mod.Coefficients{2,1}*X(:,3);
plot(X(:,3), y, 'r', 'LineWidth', 2)
ylim([0 700])
text(4, 620, ['rho = ',num2str(round(rCovar(3,5),3))])
text(4, 580, ['p = ', num2str(round(pCovar(3,5),3))])
xlabel('age')
xticks(1:10)
xticklabels(unique(score{:,4}))
ylabel('weekdaySd')

sgtitle('significant correlations between biomarkers and covariates')
clear y

