clear all
close all

displayImgs = true;
saveImgs = false;
skipLow = 2;
skipHigh = 2;
        

D = 'Samples';
% S = dir(fullfile(D,'*.jpg')); % pattern to match filenames.
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
        viscircles(centers, radii,'EdgeColor','r');
    end
    
    %%
    th = linspace(0,2*pi) ;
    for ii = 1:numel(radii)
        
        %Define the range of increasing radius sizes, from start to finish
        % - skip
        denomRange = skipLow:(radii(ii))-skipHigh; 
        
        %Genereate all y coordinates of various radius circles
        ycs = uint32(centers(ii,1)+(denomRange').*cos(th));
        %Genereate all x coordinates of various radius circles
        xcs = uint32(centers(ii,2)+(denomRange').*sin(th));
        
        %Find ys that are outside bounds
        xcsBadIdx = sum(xcs' >= nx)' ~= 0;
        xcs(xcsBadIdx,:) = [];
        ycs(xcsBadIdx,:) = [];
        %Find xs that are outside bounds
        ycsBadIdx = sum(ycs' >= ny)' ~= 0;
        xcs(ycsBadIdx,:) = [];
        ycs(ycsBadIdx,:) = [];
        
        %Caluclate (linear) index into the image for each pixel pair
        %generated above. 
        circlePixelIdx = sub2ind(size(gray_image), xcs,ycs);
        
        pixelValuesAtRs = int16(gray_image(circlePixelIdx));
        %Offset it so waveform is bipolar
        pixelValuesAtRs = pixelValuesAtRs - int16(mean(pixelValuesAtRs')');
        
        %Count the number of zero corssings.
        aboveZero = pixelValuesAtRs > 0;
        %Doing an xor like this will give a 1 when a transition of the sign
        %occurs
        diffs = xor(aboveZero(:,1:end-1),aboveZero(:,2:end));
        transitionTable = sum(diffs')';
        
        matching = transitionTable >= 3.9 & transitionTable <= 4.1;
        probability = sum(matching)/numel(transitionTable);
        
        
        strlabel = sprintf('CandidiateNumber=%d\nProbability=%3.2f\nRadius=%2.2f\n\n',...
            ii,probability*100,radii(ii));
        
        if probability > 0.8
            disp(strlabel)
        end
        if probability > 0.8 && displayImgs
            [X,Y] = meshgrid(1:ny,1:nx) ;
            % Keep only points lying inside circle
            yc = centers(ii,1)+radii(ii)*cos(th) ;
            xc = centers(ii,2)+radii(ii)*sin(th) ;
            idx = inpolygon(X(:),Y(:),yc',xc) ;
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
            if isvalid(w1h)
                viscircles(w1h, centers(ii,:), radii(ii,:),'EdgeColor','b');
            end
            if saveImgs
                saveas(iifig, strcat('Results/R_',fname,'_Circle',num2str(ii),'.jpg'));
            end
        end
    end
    if displayImgs && saveImgs
        saveas(w1h, strcat('Results/R_',fname,'_A.jpg'));
    end
end
