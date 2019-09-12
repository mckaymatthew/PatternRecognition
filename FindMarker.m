clear all
close all

displayImgs = false;
saveImgs = false;

D = 'Samples';
S = dir(fullfile(D,'Visible_IR_Daylight1PM.jpg')); % pattern to match filenames.
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

            pixelValuesAtRs = zeros(int32(radii(ii))-5-5,numel(th));
            transitionTable = zeros(int32(radii(ii))-5-5,1);
            for rDenom = 5:int32(radii(ii))-5
                radius = radii(ii);
                radiusNew = radius-double(rDenom);
                xc = uint32(centers(ii,1)+(radiusNew)*cos(th)) ; 
                yc = uint32(centers(ii,2)+(radiusNew)*sin(th)) ; 
                pixelValueAtR = zeros(numel(xc),1);
                for i = 1:numel(xc)
                    %If circle goes off image edge, handle it.
                    if yc(i) <= nx
                        c = gray_image( yc(i), xc(i));
                        pixelValueAtR(i) = c; 
                    end
                end
                pixelValueAtR = pixelValueAtR - mean(pixelValueAtR); 
                pixelValuesAtRs(rDenom+1-5,:) = pixelValueAtR;
                transitionTable(rDenom+1-5) = length(find(diff(pixelValueAtR > 0)));
            end
            
            
            matching = transitionTable == 4;
            nonMatching = transitionTable ~= 4;
            probability = sum(matching)/numel(transitionTable);
            
            
            strlabel = sprintf('CandidiateNumber=%d\nProbability=%3.2f\nRadius=%2.2f\n\n',...
                ii,probability*100,radii(ii));
            
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
