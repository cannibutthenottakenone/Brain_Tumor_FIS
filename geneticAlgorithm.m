function [weights,thresh, generations] = geneticAlgorithm(fis, images, populationN, minImprovement, eliteN, mutRate)
%GENETICALGORITHM Summary of this function goes here

% initialization
rules=length(fis.rules);

for i=1:populationN
    members(i)=PopulationMember(rand(1,rules), rand);
end

lastHiScore=Inf;
generations=0;

while true
    generations=generations+1;
   
    % evaluation
    for i=1:populationN
        members(i).evaluate(fis, images)
    end

    generationScore=sum([members.errors]);
    generationalImprovement=lastHiScore-generationScore;
    if generationalImprovement<minImprovement
        break
    end
    lastHiScore=generationScore;

    %sort based on errors (lower better)
    [~, idx] = sort([members.errors]);
    members = members(idx);

    %save elite
    if eliteN>0
        elites=members(1:eliteN);
    else
        elites=[];
    end

    % variation
    children=[];
    %recombination
    for i=1:floor((populationN-eliteN)/2)
        newWeight=zeros(2,rules);
        newThresh=zeros(2,1);

        picked=randi(populationN,2,1);
        [newWeight(1,:),newWeight(2,:)]=wholeArithmeticRecombination(members(picked(1)).weights,members(picked(2)).weights,rand);
        [newThresh(1),newThresh(2)]=wholeArithmeticRecombination(members(picked(1)).thresh,members(picked(2)).thresh,rand);

        children=[children; PopulationMember(newWeight(1,:),newThresh(1)); PopulationMember(newWeight(2,:),newThresh(2))];
    end

    %mutation
    for i=1:length(children)
        if rand<mutRate
            children(i).mutate()
        end
    end

    %replacement
    members=[elites;children];

end

% find the best member and return it's parameters
[~, idx] = min([members.errors]);
best = members(idx);

weights=best.weights;
thresh=best.thresh;

end

