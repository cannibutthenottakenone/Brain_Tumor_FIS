classdef MRImage < handle
    %MRIMAGE Class for MR image data for tumor segmentation project

    properties
        % basic properties
        image % basic data, dxd bw matrix
        label % solution to the problem

        d % dimension
        max
        min

        %-------
        % information for segmentation properties
        medianValue
        slidingMean
        slidingKurt
        gradientMods
        distanceFMedian

        floodedMap
    end

    methods
        function obj = MRImage(image)
            %MRIMAGE Construct an instance of this class
            obj.image = image;
            obj.d=size(image);
            obj.max=max(image(:));
            obj.min=min(image(:));
        end

        function obj=preprocessing(obj, windowSizes)
            %Preprocessing: populates properties reloative to the image
            %information before feeding it in the FIS

            linData=obj.image(:);
            linData=linData(linData>0); %removing background

            % - median data
            obj.medianValue=median(linData);

            % - distance from median
            obj.distanceFMedian=rescale(obj.image/obj.medianValue);

            % - mean in sliding window
            obj.slidingMean=movmean(obj.image, windowSizes.mean, "Endpoints",0);

            % - gradient in sliding window on a modified version without background
            % ('crisper' than variance)
            unblackData=obj.image;
            unblackData(unblackData==0)=obj.medianValue; %to remove strong values around the brain
            [FX,FY] = gradient(unblackData);
            modules = sqrt((FX.^2)+(FY.^2));
            obj.gradientMods = rescale(modules);

            % - kurtosis in sliding window
            obj.slidingKurt=rescale(movKurtosis(unblackData, windowSizes.kurt, obj.medianValue));
        end

        function seedCoordinates=seedAndFlood(obj, settings)
            % detect if there is at least one relevant point intensity wise
            if obj.max/obj.medianValue<settings.recognizeSeed
                obj.floodedMap=zeros(obj.d);
                seedCoordinates=[1,1];
                return
            end
            
            % select a point (more or less in the center)
            seedMap=single(obj.image==obj.max); %conversion from logical to single floating number of the map
            [seedsy, seedsx]=collapseCoordinates(seedMap);
            seedCoordinates=round([median(seedsy), median(seedsx)]);

            % identify the borders
            borders=obj.gradientMods>settings.borderThresh;
            background=obj.image==0;
            borders=borders+background; % adds the background as a wall to avoid recursion problems

            % flood
            obj.floodedMap=floodFromSeed(borders, seedCoordinates, settings.wallDistanceThresh);
            obj.floodedMap=obj.floodedMap==1; %because the flooded map contains 1 for visited and flooded places and -1 for visited but not flooded places

            %check if the flooding result is valid
            
            % if the number of elements flooded is greater than 50% of the
            % brain (obtained as number of cells minus background cells)
            % then the flood leaked everywhere
            if sum(obj.floodedMap(:))>0.5*(obj.d(1)*obj.d(2) - sum(background(:)))
                obj.floodedMap=[];
            end

            %accounting for the walldistance could be added (expand in each
            %direction by wallDistanceThresh)
        end
    end
end

