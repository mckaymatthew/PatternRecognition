clear all
close all

displayImgs = true;
saveImgs = false;

D = 'Samples';
S = dir(fullfile(D,'VisibleIRDaylight1PM.jpg')); % pattern to match filenames.
for k = 1:numel(S)
    [path, fname] = fileparts(S(k).name);
    F = fullfile(D,S(k).name);
    img = imread(F);
    
    [nx,ny,d] = size(img) ;
    img = imcrop(img, [0 140 ny nx-300]);
    [nx,ny,d] = size(img);
    gray_image = rgb2gray(img);
    
    W1 = gradientweight(gray_image, 1.5, 'RolloffFactor', 2, 'WeightCutoff', 0.25);
    [centers, radii, metric] = imfindcircles(W1,[18 54],'ObjectPolarity','dark','method','TwoStage');
    
    if displayImgs
        figure;
        imshow(W1);
        title(S(k).name);
        w1h = gca;
        viscircles(centers, radii*1.2,'EdgeColor','r');
    end
    
    
    [X,Y] = meshgrid(1:ny,1:nx) ;
    th = linspace(0,2*pi) ;
    
    %%
    for ii = 1:numel(radii)
        if ii ~= 1
            continue;
        end
        
        skipLow = 2;
        skipHigh = 1;
        pixelValuesAtRs = zeros(int32(radii(ii))-skipLow-skipHigh,numel(th));
        transitionTable = zeros(int32(radii(ii))-skipLow-skipHigh,1);
        denomRange = skipLow:(radii(ii))-skipHigh;
        for rDenom = denomRange
            radiusNew = radii(ii)-double(rDenom);
            xc = uint32(centers(ii,1)+(radiusNew)*cos(th)) ;
            yc = uint32(centers(ii,2)+(radiusNew)*sin(th)) ;
            
            pixelValueAtR = int16(gray_image(sub2ind(size(gray_image), yc,xc)));
            pixelValueAtR = pixelValueAtR - mean(pixelValueAtR);
            pixelValuesAtRs(rDenom+1-skipLow,:) = pixelValueAtR;
            transitionTable(rDenom+1-skipLow) = length(find(diff(pixelValueAtR > 0)));
        end
%         pixelValuesAtRs = (pixelValuesAtRs' - mean(pixelValuesAtRs'))';
        
        
        matching = transitionTable >= 3.9 & transitionTable <= 4.1;
        probability = sum(matching)/numel(transitionTable);
        
        
        strlabel = sprintf('CandidiateNumber=%d\nProbability=%3.2f\nRadius=%2.2f\n\n',...
            ii,probability*100,radii(ii));
        
        %             probability = 1
        if probability > 0.8
            disp(strlabel)
        end
        if probability > 0.8 && displayImgs
            % Keep only points lying inside circle
            xc = centers(ii,1)+radii(ii)*cos(th) ;
            yc = centers(ii,2)+radii(ii)*sin(th) ;
            idx = inpolygon(X(:),Y(:),xc',yc) ;
            bound = [centers(ii,1)-radii(ii) centers(ii,2)-radii(ii) radii(ii)*2 radii(ii)*2];
            
            for i = 1:d
                I1 = img(:,:,i) ;
                I1(~idx) = 255 ;
                imgc(:,:,i) = I1 ;
            end
            
            wImage = W1;
            wImage(~idx) = 255;
            gImage = gray_image;
            gImage(~idx) = 255;
            
            
            cImage =  imcrop(imgc,bound);
            wImage =  imcrop(wImage,bound);
            gImage =  imcrop(gImage,bound);
            gImageEdge = edge(gImage,'canny');
            plotrows = 5;
            figure;
            subplot(plotrows,2,1);
            imshow(wImage);
            subplot(plotrows,2,2);
            imshow(edge(wImage,'canny'));
            subplot(plotrows,2,3);
            imshow(cImage);
            subplot(plotrows,2,5);
            imshow(gImage);
            subplot(plotrows,2,6);
            imshow(gImageEdge);
            subplot(plotrows,2,[9 10]);
            ylabel('Intensity');
            xlabel('Pixel');
            plot(pixelValuesAtRs');
            subplot(plotrows,2,[7 8]);
            
            plot(transitionTable);
            ylabel('Zero Crossings');
            xlabel('Radius (Reversed)');
            
            annotation('textbox', [0, 0.5, 0, 0], 'string',strlabel);
            iifig = gca;
            viscircles(w1h, centers(ii,:), radii(ii,:),'EdgeColor','b');
            if saveImgs
                saveas(iifig, strcat('Results/R_',fname,'_Circle',num2str(ii),'.jpg'));
            end
        end
    end
    if displayImgs && saveImgs
        saveas(w1h, strcat('Results/R_',fname,'_A.jpg'));
    end
end
