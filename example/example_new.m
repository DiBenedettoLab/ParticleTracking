%% Particle tracking code example
% this example is for an .avi file but can be changed for tiff stacks or
% .seq files easily

close all
clear
clc

%% Filename
addpath(genpath('../src/'))
inputname = 'data/testtracks.avi';

%% Test image (not necessary)
vid = VideoReader(inputname);

%read one frame
frame_num = 1;
img = read(vid,frame_num);

% plot one frame
figure
imagesc(img)

% use datatips to click around on the image to get an idea of the size of
% the particle and its intensity difference with the background (will be
% helpful for the next part)

%% Particle finding

% INPUTS: (you may need to adjust these)

% threshold (difference in intensity between particle and background)
threshold = 40;

%background (name of background image if it exists. here it is just empty)
bground_name = [];

%minimum area (minimum area of particle to look for in pixels.) for large
%particles, this should be greater than 1.
minarea = 2;

%invert (=0 if light particles on dark surface. =1 if dark particles on
%light surface)
invert = 0;

%image frames to look at. if left empty [], then does all frames
%vector of length two with first and last index designated
framerange=[];%[ 11 20];

[centroids] = Particle_Finder(inputname,threshold,framerange,[],bground_name,minarea,invert,0);

%% Particle tracking

%maximum displacement (max allowable pixels a particle is allowed to move
%between frames.) (this is especially important when there are many
%particles in the field of view)
max_disp = 8;

[vtracks,ntracks,meanlength,rmslength] = Predictive_Tracker(centroids,'max_disp',max_disp);


%% Another useful function
% this function converts the vtracks into vectors (helpful for when there
% are many particles). 
[u,v,x,y,t]=Velocities(vtracks);

%% Plot outputs
% just plot particles
figure(2)

if isempty(framerange)
images2plot=1:vid.NumFrames; %plots all frames
else
images2plot=framerange(1):framerange(2);
end

for i=images2plot
    
    hold off
    imagesc(read(vid,i))
    hold all
    
    part = find(i==t);
    scatter(x(part),y(part),'ro')
    
    pause(.1) %seconds to pause
    
end


%% Plot tracks on images
plot_tracks_avi(inputname,vtracks,framerange)


%% Clean tracks

% remove tracks that are short in time and space

% Need to change the following 3 parameters to adjust cleaning:
min_length = 3; %frames (could set to 0 to keep all)
min_y_dist = 0; %pixels (min distance moved in y)
min_x_dist = 0; %pixels (min distance moved in x)

good = ones(length(vtracks),1);

for i=1:length(vtracks)
    
    % check if its long enough in time
    t_dist = length(vtracks(i).T)< min_length;
    
    % check if its long enough in space
    x_dist = ( max(vtracks(i).X)-min(vtracks(i).X) ) < min_x_dist;
    y_dist = ( max(vtracks(i).Y)-min(vtracks(i).Y) ) < min_y_dist;
    
    % if t_dist, x_dist, or y_dist is "true", then throw out the track
    
    if t_dist || x_dist || y_dist
        good(i) = 0;
    end
    
end

vtracks_cleaned=vtracks(logical(good));

disp(['Total tracks before cleaning: ' num2str(numel(vtracks),'%.0f')])
disp(['Total tracks too short: ' num2str(sum(~good),'%.0f')])
disp(['Total tracks after cleaning: ' num2str(sum(good),'%.0f')])



%% saving info

% will save your vtracks variable to a file
% make sure to either rename it or move it so you don't overwrite a
% previous version

save('testtracks_saved.mat','vtracks_cleaned')
