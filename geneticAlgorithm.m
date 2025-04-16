function [weights,thresh, generations] = geneticAlgorithm(fis, images, populationN, minImprovement, eliteN, mutRate)
%GENETICALGORITHM Summary of this function goes here

% initialization
rules=length(fis.rules);

for i=1:populationN
    members(i)=PopulationMember(rand(1,rules), rand);
end

lastHiScore=Inf;
improvements=[];
bestFit=[]; %should go down
generations=0;

oldM=[]; % saves old generation in case of next generation is much worse

while true
    generations=generations+1;
    fprintf("Generation:%d\n\n", generations)
   
    % evaluation
    for i=1:length(members) %more than population because after recombination we get more, half of which will then be killed
        fprintf("Evaluation of genome:%d\n", i)
        members(i).evaluate(fis, images)
    end

    generationScore=sum([members.errors]);
    generationalImprovement=lastHiScore-generationScore;
    improvements(generations)=generationalImprovement;
    if generationalImprovement<-500*populationN % new generation is much worse
        generationalImprovement=0
        improvements(generations)=0;
        members=oldM;
        fprintf("New generation was so bad i threw it away\n")
    end
    bestFit(generations)=min([members.errors]);
    fprintf("improvement wrt last generation:%d less errors\n\n", generationalImprovement)
    if sum(improvements(1:min(length(improvements),10))<minImprovement*populationN)>=5
        break
    end
    lastHiScore=generationScore;

    %sort based on errors (lower better)
    [~, idx] = sort([members.errors]);
    members = members(idx);

    % save all members in case of next generation is much worse
    oldM=members;

    %parent selection, the best N are saved
    members=members(1:populationN);

    %save elite
    if eliteN>0
        elites=members(1:eliteN);
    else
        elites=[];
    end

    % variation
    children=[];
    %recombination
    for i=1:populationN-eliteN
        newWeight=zeros(2,rules);
        newThresh=zeros(2,1);
        
        %sort of random tournament, pick random 4, recombine the 2 best
        candidates=sort(randi(populationN,4,1));
        picked=candidates([1,2]);

        [newWeight(1,:),newWeight(2,:)]=wholeArithmeticRecombination(members(picked(1)).weights,members(picked(2)).weights,rand);
        [newThresh(1),newThresh(2)]=wholeArithmeticRecombination(members(picked(1)).thresh,members(picked(2)).thresh,rand);

        children=[children, PopulationMember(newWeight(1,:),newThresh(1)), PopulationMember(newWeight(2,:),newThresh(2))];
    end

    %mutation
    for i=1:length(children)
        if rand<mutRate
            children(i).mutate()
        end
    end

    %replacement
    members=[elites,children];

end

% find the best member and return it's parameters
[~, idx] = min([members.errors]);
best = members(idx);

weights=best.weights;
thresh=best.thresh;

end

