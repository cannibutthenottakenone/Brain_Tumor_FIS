classdef PopulationMember < handle
    %POPULATIONMEMBER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        weights
        thresh
        errors
    end

    methods
        function obj = PopulationMember(weights,thresh)
            obj.weights=weights;
            obj.thresh=thresh;
            obj.errors=0;
        end

        function evaluate(obj,fis, images)
            % apply weights to fis
            for i=1:length(obj.weights)
                fis.rules(i).Weight=obj.weights(i);
            end

            for i=1:length(images)
                % computation of inferences
                if isempty(images(i).floodedMap)
                    flood=ones(images(i).d)*-1;
                else
                    flood=images(i).floodedMap;
                end
                input=[images(i).slidingMean(:),images(i).gradientMods(:),flood(:),images(i).slidingKurt(:),images(i).distanceFMedian(:)];

                result=reshape(evalfis(fis, input), images(i).d);
                result=result>obj.thresh;

                % computation of scores
                matrixScore=images(i).label-result;
                obj.errors=obj.errors+nnz(matrixScore(:));
            end
        end

        function mutate(obj) %implements non-uniform mutation
             obj.weights=obj.weights+normrnd(0,0.03,1,length(obj.weights));

             obj.weights(obj.weights>1)=1;
             obj.weights(obj.weights<0)=0;

             obj.thresh=obj.thresh+normrnd(0,0.03);
        end
    end
end

