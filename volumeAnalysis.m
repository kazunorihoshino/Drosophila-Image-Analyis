% volumeAnalysis 07252024


%% Step 0: file information

clc       % close all
clear     % clear all
close all

%file names
imgpath   = '/Users/stellacho/Desktop/OA-/';
datapath  = '/Users/stellacho/Desktop/OA-/datatest/';
videopath = '/Users/stellacho/Desktop/OA-//movtest/';

filename = 'Control1.tif';


%layer information
zlayers = 55;     %number of zslice
Nf = 60; %number of steps

samplelayer = 26;  % Choose a decently focused layer

% size of window
halfwidth = 192;


%% Step 1: Select a sample

mkdir(datapath);
mkdir(videopath);

filenamefull=strcat(imgpath,filename);

info = imfinfo(filenamefull);
totalNf = numel(info)


% bottom five and top four layers will be wiped  for margins %

imgs = zeros (halfwidth *2, halfwidth *2, zlayers, Nf, 'uint8');
testimg = imread(filenamefull,samplelayer);

% you can skip this if you did this before

button = questdlg('have you already located the sample? if this is the first time opening the saemple, click No','find center','Y','N','Y')

if button == 'Y'
    locationfilename=strcat(datapath,filename,'.location.mat');
    load(locationfilename);
    angle = location(1);
    m_y  = floor(location(2));
    m_x  = floor(location(3));
    
    rotate = imrotate(testimg,-angle);

    fullwidth = halfwidth * 2 + length(rotate);
    largeimg = zeros(fullwidth,fullwidth,'uint8');

    largeimg(halfwidth + 1: halfwidth + length(rotate),halfwidth + 1: halfwidth + length(rotate)) = rotate;
    
    crop = largeimg(m_y - halfwidth:m_y + halfwidth-1,m_x - halfwidth : m_x + halfwidth-1);
    
    if location(4) == -1
        crop = imrotate(crop,180);
    end
    
    figure;imshow(crop) 
    
else
    
    waitfor(msgbox('Click the top and the bottom ends.'))
    [P_Y,P_X] = clickpoints(testimg,2);

    angle = (180/pi) * atan( (P_Y(2)-P_Y(1)) / (P_X(2)-P_X(1)) );
    location(1) = angle;
    rotate = imrotate(testimg,-angle);
   
    imshow(rotate)
    
    
    fullwidth = halfwidth * 2 + length(rotate);
    largeimg = zeros(fullwidth,fullwidth,'uint8');

    largeimg(halfwidth + 1: halfwidth + length(rotate),halfwidth + 1: halfwidth + length(rotate)) = rotate;
    waitfor(msgbox('Click the two points again. Top first.'))
    [P_Y,P_X] = clickpoints(largeimg,2);

    m_y = floor(mean(P_X)) 
    m_x = floor(mean(P_Y))
    
    location=[angle,mean(P_X),mean(P_Y)];
    

    crop = largeimg(m_y - halfwidth:m_y + halfwidth-1,m_x - halfwidth : m_x + halfwidth-1);

    
    if (P_X(1)>P_X(2))
        crop = imrotate(crop,180);
        location(4) = -1;
    else
        location(4)=1;
    end
    close all;
    figure;imshow(crop)
    
    locationfilename=strcat(datapath,filename,'.location.mat');
    save(locationfilename,'location'); 
    
end


waitfor(msgbox('Redo this step if the sample is not in the middle.'))


%% Step 2: file reading. It takes time.
       
close all;
zs = 0;

for i = 1:1:Nf
    for j = 1:1:zlayers            
        fullimage = imread(filenamefull,(zs + j) + zlayers*(i-1));
        rotate = imrotate(fullimage,-angle);
        
        largeimg = zeros(fullwidth,fullwidth,'uint8');
        largeimg(halfwidth + 1: halfwidth + length(rotate),halfwidth + 1: halfwidth + length(rotate)) = rotate;
        
        imgs(:,:,j,i) = largeimg(m_y - halfwidth:m_y + halfwidth-1,m_x - halfwidth : m_x + halfwidth-1);
            
        if location(4)==-1
            imgs(:,:,j,i) = imrotate(imgs(:,:,j,i),180);
        end
    end
    i
    figure;imshow(imgs(:,:,samplelayer,i));
end

%% Step 3: Find the rim.


nodenx = 64 + 1;   %number of nodes in X direction
nodeny = 64 + 1;   %number of nodes in Y direction
nodenz = zlayers;   %number of nodes in Z direction. We only track the surface 12/24/21

areanx = nodenx - 1;   %number of nodes in X direction
areany = nodeny - 1;   %number of nodes in Y direction
areanz = nodenz - 1;   %number of nodes in Z direction

pixnx = (halfwidth * 2)/ areanx;   %size of an element in X
pixny = (halfwidth * 2)/ areany;   %size of an element in Y
pixnz = 1;                         % 


imgwithcolorrects = zeros (halfwidth *2,halfwidth *2,3,nodenz, Nf,'uint8');

areas       = zeros(areanx, areany, areanz,'uint8');
areas_f     = zeros(areanx, areany, Nf,'uint8');
areas_Nf    = zeros(areanx, areany, areanz,'uint8');


        BW = imbinarize(imgs(:,:,samplelayer,1), 0.035); 


        figure;imshow(BW)
        
        BW2 = imfill(BW,'holes');
        
        imshow(BW2);
        
        

        for ii = 1:1:1    %number of trims
            BW3 = bwperim(BW2,8);
            BW2 = BW2 -BW3;
        end

        
        
       figure;imshow(BW2);
        
       
       
% find the layer with the most focused rim. You can skip this
close all;

button = questdlg('Find the layer with the most focused rim. If you have alreay done this, click No','find center','Y','N','Y')
   

if button == 'Y'

    for k = 1:1:areanz   

       disp = imgs(:,:,k,1); 
       figure
       imshow(disp)      
    end
    
else
    k = uint8(location(5))

end


%% Step 4: Create mask files 1
% create BW area maps
close all;


if button == 'Y'
    kanswer = inputdlg('Enter the layer number or cancel to use the previous value (if there is no previous value and cancerl, you will see an error).'); %for example 31
    if size(kanswer) ~= [0 0]
        location(5)=str2num(kanswer{1});
    end
    k = uint8(location(5))
    locationfilename=strcat(datapath,filename,'.location.mat');
    save(locationfilename,'location'); 
end

for f=1:1:Nf

   BW = imbinarize(imgs(:,:,k,f), 0.035); 


   BW2 = imfill(BW,'holes');

   for ii = 1:1:1    %number of trims
        BW3 = bwperim(BW2,8);
        BW2 = BW2 -BW3;
   end    
        
   for i = 1:1:areanx
        for j = 1:1:areany
            part = BW2(pixnx * (i-1) +1 : pixnx * i , pixny * (j-1) +1 : pixny * j);
        
            part = reshape (part, [pixnx * pixny,1]);
        
            if mean(part) > 0.6     %we put threshold at anywhere 0.4-0.6 or something... 
        
                BW2(pixnx * (i-1) +1 : pixnx * i , pixny * (j-1) +1 : pixny * j) = 1;
                areas_f(i,j,f) = 255;
            else
                BW2(pixnx * (i-1) +1 : pixnx * i , pixny * (j-1) +1 : pixny * j) = 0;
                areas_f(i,j,f) = 0;
            end
            
        end
   end  
        
    maskfilename = strcat(datapath,'BW_f_',num2str(f,'%03d'),'.bmp');    
    imwrite(areas_f(:,:,f),maskfilename,'BMP')
    
    
end

maskfilename = strcat(datapath,'BW_f_000.source.bmp');    
imwrite(areas_f(:,:,1),maskfilename,'BMP')

maskfilename = strcat(datapath,'BW_f_',num2str(Nf,'%03d'),'.source.bmp');    
imwrite(areas_f(:,:,f),maskfilename,'BMP')


waitfor(msgbox('Modify BW_f_000.source.bmp and create BW_f_000.bmp.If you already have the file, do not make a new one'))



%% Step 5: Create mask files 2
close all;


    maskfilename = strcat(datapath,'BW_f_000.bmp');
    areas_f0 = imread(maskfilename);

for f=1:1:Nf
  maskfilename = strcat(datapath,'BW_f_',num2str(f,'%03d'),'.bmp');
  areas_f(:,:,f) = imread(maskfilename);
    
end


%

areas       = zeros(areanx, areany, areanz,'uint8');

for f=1:1:Nf
    
   for i = 1:1:areanx
        for j = 1:1:areany
            areas_f(:,:,f) = areas_f(:,:,f) .* areas_f0(:,:); % areas_f0(:,:) is a mask. (.*) is a product element by element
            
        end
   end
   
    maskfilename = strcat(datapath,'BW_f_',num2str(f,'%03d'),'.bmp');   
    imwrite(areas_f(:,:,f),maskfilename,'BMP')
end


%

for f=1:1:Nf
  maskfilename = strcat(datapath,'BW_f_',num2str(f,'%03d'),'.bmp');
  areas_f(:,:,f) = imread(maskfilename);
    
end


% contrast test and create focus map

close all;
delta_x = 28;  % Size of focusing window
delta_y = 28;


half_delta_x = floor(delta_x/2);
half_delta_y = floor(delta_y/2);


contrast =  zeros(areanx,areany,zlayers -1);
focusmap =  zeros(areanx,areany,zlayers -1,Nf,'uint8');
heightmap = zeros(areanx,areany,Nf,'uint8');
topimgs   = zeros(halfwidth *2, halfwidth *2,Nf,'uint8');
focusplot = zeros(zlayers -1,1,'uint8')


for f=1:1:60
    f
    maskfilename = strcat(datapath,'BW_f_',num2str(f,'%03d'),'.bmp');
    areas_f(:,:,f) = imread(maskfilename);
        
    
    for k = 1:1:zlayers -1   
        for j = 1:1:areany
            for i = 1:1:areanx
                if areas_f(i,j,f) == 255
                    window = (imgs(1 + pixnx * (i-1)+ pixnx/2-half_delta_x : pixnx * (i-1)+ pixnx/2+half_delta_x, 1 + pixny * (j-1)+ pixny/2-half_delta_y : pixny * (j-1)+ pixny/2+half_delta_y,k,f));
                    contrast(i,j,k) =  (max(max(window))-min(min(window)));
                end
            end
        end
    end


    for j = 1:1:areany
        for i = 1:1:areanx
            if areas_f(i,j,f) == 255

              [maxcont maxii] = max(contrast(i,j,:));
              [mincont minii] = min(contrast(i,j,:));

              focusmap(i,j,:,f) = 255*(contrast(i,j,:)-mincont)/(maxcont-mincont);
              
              %focusplot(:) = focusmap(i,j,:,f) ;
              %figure;plot(focusplot)              
              
              heightmap(i,j,f) = min(50,max(10,maxii)); 
              
              topimgs(1+(i-1)*pixnx:i*pixnx,1+(j-1)*pixny:j*pixny,f) = imgs(1+(i-1)*pixnx:i*pixnx,1+(j-1)*pixny:j*pixny,heightmap(i,j,f),f);
                           
            end
         end
    end
    
    figure;imshow(topimgs(:,:,f));

    if f == 1;
        for k = 1:1:areanz   %these are elements newly added on 12/23/2021 (commented out the nodes plot above)

           %figure;imshow(imresize(focusmap(:,:,k,f),pixnx));

        end
    end
end

waitfor(msgbox('Modify BW_f_000.bmp and repeat if needed.'))


%% Step 6: Check files

    
for k = 14:1:84   
    for j = 1:1:areany
        for i = 1:1:areanx
            if heightmap(i,j,1)==k
                areas(i,j,k) = 255;
            end
        end
    end      
end


for k = 1:1:areanz                     
    maskfilename = strcat(datapath,'BW4_z',num2str(k,'%03d'),'.bmp');    
    imwrite(areas(:,:,k),maskfilename,'BMP')                  
end


%
close all;

for k=1:1:areanz
    maskfilename = strcat(datapath,'BW4_z',num2str(k,'%03d'),'.bmp');
    areas(:,:,k) = imread(maskfilename);
end


%% Step 7: Volume analysis. Count, find, and visualize the rim height
close all;

rimheight = zeros (Nf,1,'double');
rimview = zeros(areanx,areany,Nf,'uint8');

areacountN = 100;

for f = 1:1:Nf
    areacount = 0;
    for k = areanz:-1:1
        for j = 1:1:areany
            for i = 1:1:areanx
                if heightmap(i,j,f)==k
                    if areacount <areacountN
                        areacount = areacount + 1
                        rimheight(f) = rimheight(f) + k;
                        rimview(i,j,f) = 255;
                    end
                end
            end
        end
    end
    rimheight(f) = rimheight(f)/areacountN;
    figure;imshow(rimview(:,:,f));
end

% evaluate the height map
close all

volume = zeros(Nf,1,'double');
topviewarea = zeros(Nf,1,'uint32');

 maxheight = zeros(Nf,1,'uint8');
 minheight = zeros(Nf,1,'uint8');
 
 for f = 1:1:Nf
    maxheight(f) = max(max(heightmap(:,:,f)));
    minheight(f) = 50;
     for j = 1:1:areany
        for i = 1:1:areanx
            if heightmap(i,j,f)>0                
                minheight(f) = min(minheight(f),heightmap(i,j,f));  
            end
        end
     end
 end



for f = 1:1:Nf
    
    for j = 1:1:areany
        for i = 1:1:areanx
            if heightmap(i,j,f)>0
                volume(f) = volume(f) + double(rimheight(f)) - double(heightmap(i,j,f)) ;
                topviewarea(f) = topviewarea(f) + 1;
            end
        end
    end
end


figure;plot(minheight)
figure;plot(maxheight)
figure;plot(rimheight)

figure;plot(volume);


volumefilename = strcat(filename,'.volume.csv');    
writematrix(volume,volumefilename) 
