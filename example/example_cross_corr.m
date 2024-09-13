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

[u,v,x,y,t,tr,area]=Extract_vtracks(vtracks);

%% Find largest particles to test cross-corr on

minarea = 7;

big_area_tracks = find(area>minarea);

tr_big = mode(big_area_tracks);


%% get cross_corr velocity and compare to particle tracker
close all
n = vtracks(tr_big).len;
medfilt_win = 0;
windowsize = 5;

df = 1; %frames between cross-corr

yesplot = 1;

for i = 1:n-df;
fr1 = vtracks(tr_big).T(i);
fr2 = vtracks(tr_big).T(i+df);

X = vtracks(tr_big).X([i i+df]);
Y = vtracks(tr_big).Y([i i+df]);

U = vtracks(tr_big).U([i i+df]);
V = vtracks(tr_big).V([i i+df]);


im1 = read(vid,fr1);
im2 = read(vid,fr2);


[u_corr(i,:),u_guess(i,:),corr_peaks(i)] = Particle_Velocity_Cross_Corr(X,Y,im1,im2,windowsize,medfilt_win,yesplot);

u_tracker(i,:) = [mean(U) mean(V)];

pause(.1)

end
u_corr = u_corr/df;
u_guess = u_guess/df;


figure
plot(u_corr)
hold all
plot(u_tracker,'--')

pause
%% increase df to minimize peaklocking
figure
n = vtracks(tr_big).len;
medfilt_win = 0;
windowsize = 5;

df = 2; %frames between cross-corr

yesplot = 0;

for i = 1:n-df;
fr1 = vtracks(tr_big).T(i);
fr2 = vtracks(tr_big).T(i+df);

X = vtracks(tr_big).X([i i+df]);
Y = vtracks(tr_big).Y([i i+df]);

U = vtracks(tr_big).U([i i+df]);
V = vtracks(tr_big).V([i i+df]);


im1 = read(vid,fr1);
im2 = read(vid,fr2);


[u_corr(i,:),u_guess(i,:),corr_peaks(i)] = Particle_Velocity_Cross_Corr(X,Y,im1,im2,windowsize,medfilt_win,yesplot);

u_tracker(i,:) = [mean(U) mean(V)];

pause(.1)

end
u_corr = u_corr/df;
u_guess = u_guess/df;


figure
plot(u_corr)
hold all
plot(u_tracker,'--')