clear all
format long
R=VideoReader ('visiontraffic.avi'); %creates object v to read video data from the file named filename
X_Data=Matrix_video(R,171,200);%Put the 30 frames in a matrix
[n,m]=size(X_Data);
X = double(X_Data); % matricea cu cele 30 de frame-uri din video

L1 = rand(n, m); 
L2 = randn(n, m); 
S1 = rand(n, m); 
S2 = rand(n, m); 
lagr = 0.0002 * eye(n, m); 

p = 0.007;
lambda = 0.0002; 
epsi = 0.1;

reziduuri = [];
iter = 0;

while rank(L2) > 1 && iter < 200
    reziduuri = [reziduuri, norm(L2 + S2 - X)];
    norm(L2 + S2 - X)

    L1 = L2; %am copiat valorile vechi pentru pasul urmator
    S1 = S2;
    Y = X - S1 - (1 / p) * lagr;
   
    Y = X - S1 - (1 / p) * lagr; %actualizez L
    L2 = SVT(Y, 1/p); 

    Z = X - L2 - (1 / p) * lagr; %actualizez S
    S2 = sign(Z) .* max(abs(Z) - lambda / p, 0); 

    lagr = lagr + p * (L2 + S2 - X);
    iter = iter + 1;
end
Lnew = L2;
Snew = S2;
semilogy(reziduuri);

%% Decomentati pentru a crea un video, dupa ce ati implementat algoritmul de optimizare
% %% Loop through the image sequence and write each image as a frame to the video
% % Specify the output video file name and format
outputVideo = VideoWriter('rezultat.avi');
outputVideo.FrameRate = 10; % You can set your desired frame rate here

% Open the video file for writing
open(outputVideo);


for j = 1:m
 Img_bg=Lnew(:,j); % background image
 Img = Snew(:,j);  % foreground image
 Img_org= X_Data(:,j);
    writeVideo(outputVideo, [reshape(Img_org,360,640), uint8(reshape(Img_bg,360,640)), reshape(Img,360,640)]);
end

% Close the video file
close(outputVideo);