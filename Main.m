% I suggest using the live script version of this file (LiveMain.mlx) to have a
% better visualization of the steps of the program

%% settings

windowSizes=struct("mean",4,"kurt",4); %sizes of sliding windows
floodingSettings=struct("borderThresh",0.25,"wallDistanceThresh", 2, "recognizeSeed",1.3); %settings for flooding:
% - borderThresh: threshold above which a normalized gradient is considered a
%border
% - wallDistanceThresh: distance from a wall that will stop the flooding
% - to be able to find an acceptable seed the max intensity must be at
% least {recognize seed} times the median intensity

%% Data files identification 
imageDir = fullfile(pwd, "BraTS","sampleBraTSTestSetValid", "imagesTest");
labelDir = fullfile(pwd, "BraTS","sampleBraTSTestSetValid", "labelsTest");
imageFiles = dir(fullfile(imageDir, "*.mat"));
labelFiles = dir(fullfile(labelDir, "*.mat"));

% Extract file names
imageSubjects = string(erase({imageFiles.name},".mat"));
labelSubjects = string(erase({labelFiles.name},".mat"));

% Visualize available files
fprintf("Soggetti con immagini disponibili: %s\nSoggetti con etichette disponibili: %s\n", strjoin(imageSubjects, ", "), strjoin(labelSubjects, ", "))
clear imageFiles labelFiles

%% slices extraction
slicesNumber=10;
rng(0,'twister');

slices=strings(slicesNumber,2);

slices(:,1)=imageSubjects(randi([1,length(imageSubjects)], slicesNumber, 1));
slices(:,2)=randi([40,112],slicesNumber, 1);

%% Data Loading
images=[];

for i=1:slicesNumber
    % Data loading (could be optimized since it loads the same file so many times)
    dataFileName = fullfile(imageDir,slices(i,1));
    labelFileName = fullfile(labelDir,slices(i,1));
    data = load(dataFileName);
    label = load(labelFileName);
    data = data.cropVol(:,:,str2double(slices(i,2)),1);
    label = label.cropLabel(:,:,str2double(slices(i,2)));

    image=MRImage(data);
    image.label=label;

    %object creation and assignment to matrix
    images=[images,image];
end

imagesL=length(images);


%loading of test slides
tslices=strings(10,2);
tslices(:,1)=imageSubjects(randi([1,length(imageSubjects)], 10, 1));
tslices(:,2)=randi([40,112],10, 1);

for i=1:10
    dataFileName = fullfile(imageDir,tslices(i,1));
    labelFileName = fullfile(labelDir,tslices(i,1));
    data = load(dataFileName);
    label = load(labelFileName);
    data = data.cropVol(:,:,str2double(tslices(i,2)),1);
    label = label.cropLabel(:,:,str2double(tslices(i,2)));

    image=MRImage(data);
    image.label=label;

    %object creation and assignment to matrix
    tImages(i)=image;
end

%cleanup
clear data dataFileName i image imageDir imageSubjects label labelDir labelFileName labelSubjects slicesNumber

%% Compute various information on the image
fprintf("preprocessing...\n")
for i=1:imagesL
    images(i).preprocessing(windowSizes); %computes all values displayed underneath
end
clear i ans

%% Seeding and flooding to identify a cancerous area based on brightest point seed and gradient walls
coordinates=zeros(imagesL,2);
for i=1:imagesL
    coordinates(i,:)=images(i).seedAndFlood(floodingSettings);
end

%% Application to the FIS and genetic alorithm to select the best weights
fis=readfis("SugenoImageRecognition.fis");
%% Evaluating images before weight tuning
results=cell(imagesL, 1);
for i=1:imagesL
    if isempty(images(i).floodedMap)
        flood=ones(images(i).d)*-1;
    else
        flood=images(i).floodedMap;
    end
    input=[images(i).slidingMean(:),images(i).gradientMods(:),flood(:),images(i).slidingKurt(:),images(i).distanceFMedian(:)];
    result=reshape(evalfis(fis, input), images(i).d);
    results{i}=result; %should be split to binary at 0.5, but not done for visualization purposes 
end
clear flood i input result

%% Tuning
[bestWeights, bestThresh, genNumber]=geneticAlgorithm(fis, images, 5, 300, 3, 0.9); %it takes a while
fprintf("best weights: %s.\nbest thershold:%f.",strjoin(string(bestWeights), ", "), bestThresh)

%% Final FIS
%replacing weights
for i=1:length(bestWeights)
    fis.rules(i).Weight=bestWeights(i);
end

%computing inferences
resultsf=cell(imagesL, 1);
for i=1:imagesL
    if isempty(images(i).floodedMap)
        flood=ones(images(i).d)*-1;
    else
        flood=images(i).floodedMap;
    end
    input=[images(i).slidingMean(:),images(i).gradientMods(:),flood(:),images(i).slidingKurt(:),images(i).distanceFMedian(:)];
    result=reshape(evalfis(fis, input), images(i).d);
    resultsf{i}=result>bestThresh;
end

%% Evaluation
%preprocessing
tImagesL=length(tImages);

for i=1:tImagesL
    tImages(i).preprocessing(windowSizes);
    tImages(i).seedAndFlood(floodingSettings);
end
%computing inferences
resultTest=cell(tImagesL, 1);
errors=zeros(tImagesL, 1);
for i=1:tImagesL
    if isempty(tImages(i).floodedMap)
        flood=ones(tImages(i).d)*-1;
    else
        flood=tImages(i).floodedMap;
    end
    input=[tImages(i).slidingMean(:),tImages(i).gradientMods(:),flood(:),tImages(i).slidingKurt(:),tImages(i).distanceFMedian(:)];
    result=reshape(evalfis(fis, input), tImages(i).d);
    resultTest{i}=result>bestThresh;
    errorMatrix=tImages(i).label-resultTest{i};
    errors(i)=nnz(errorMatrix);
end
fprintf("average errors per picture:%f pixels",sum(errors)/tImagesL)

figure
for i=1:tImagesL
    subplot(2,tImagesL,i)
    imshow(tImages(i).label)
    title(sprintf("label %d",i))
    subplot(2,tImagesL,tImagesL+i)
    imshow(resultTest{i})
    title(sprintf("predition %d",i))
end