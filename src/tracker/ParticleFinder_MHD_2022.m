function [x,y,t,ang] = ParticleFinder_MHD(inputnames,threshold,framerange,outputname,bground_name,arealim,invert,noisy)
% Usage: [x,y,t,ang] = ParticleFinder(inputnames,threshold,[framerange],[outputname],[bground_name],[arealim],[invert],[noisy])
% Given a movie of particle motions, ParticleFinder identifies the
% particles, returning their positions, times, and orientations in x, y,
% and t, respectively. The movie must be saved as a series of image files,
% an image stack in .tif or .gif format, or an uncompressed .avi file;
% specify the movie in "inputnames" (e.g., '0*.png' or 'stack.tif', or
% 'movie.avi'). To be identified as a particle, a part of the image must
% have brightness that differs from the background by at least "threshold".
% If invert==0, ParticleFinder seeks particles brighter than the
% background; if invert==1, ParticleFinder seeks particles darker than the
% background; and if invert==-1, ParticleFinder seeks any sort of contrast.
% The background is read from the file "bground_name"; see BackgroundImage.
% Frames outside the range specified by the two-element vector "framerange"
% are ignored. If arealim==1, ParticleFinder seeks single-pixel particles
% by comparing brightness to adjacent pixels (fast and good for small
% particles); otherwise ParticleFinder seeks particles having areas bounded
% by the two elements of the vector "arealim" (in square pixels; this
% method is better for tracking large particles). If "outputname" is not
% empty, particle positions are also saved as a binary file of that name.[]

% The file begins with a/Users/mdiben/Desktop long int giving the frame count, then each frame
% begins with a long int giving its particle count, and continues with
% double-precision floats giving the x and y coordinates of each particle.
% If noisy~=0, the movie is repeated with particle locations overlaid. If
% noisy>1, each movie frame is also saved to disk as an image. See also
% BackgroundImage.m and PredictiveTracker.m. Requires
% read_uncompressed_avi.m for use with .avi movies.

% Written 20 October 2011 by Doug Kelley, largely based on
% PredictiveTracker.m.
% Renamed ParticleFinder and incorporated FindParticles function 27 October
% 2011.
% Fixed plotting bug in bug plotting 15 November 2011.
% Added "ang" output (particle orientation) 18 November 2011.
% Added "framerange" input 28 November 2011.
% Updated to use weighted centroid 2 December 2011.
% Included invert==-1 option 22 February 2012.
% Made compatible with tiff & gif stacks and squelched regionprops
% divide-by-zero warning 7 March 2012.
% Thanks go to Max Gould for spotting a bug in the main loop; fixed 26
% March 2012.
% Changed 29 March 2012 to plot bars if minarea>1.
% Updated 30 March 2012 to use *weighted* orientaion and pre-allocate
% arrays.
% Updated 29 May 2012 with "arealim" input instead of "minearea".

% Next: write angle to output file!

% -=- Set defaults -=-----------------------------------------------------
framerange_default = [1 inf]; % by default, all frames
bground_name_default = 'background.tif';
noisy_default=0; % don't plot unless requested
arealim_default=1;
invert_default=-1; % by default, use absolute contrast
pausetime=1/30; % seconds to pause between frames when plotting
savedirname='particlesmovie';
barsize=8; % if interested in angle, plot particles as rods 5 px long

if nargin<2
    error(['Usage: [x,y,t,ang] = ' mfilename ...
        '(inputnames,threshold,[framerange],[outputname],' ...
        '[bground_name],[arealim],[invert],[noisy])'])
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerange_default;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('bground_name','var') || isempty(bground_name)
    bground_name=bground_name_default;
end
if ~exist('arealim','var') || isempty(arealim)
    arealim=arealim_default;
end
if ~exist('invert','var') || isempty(invert)
    invert=invert_default;
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisy_default;
end
if ~exist('outputname','var') || isempty(outputname)
    writefile=false;
else
    writefile=true;
end
%if ~exist('index','var') || isempty(index)  %%MHD
    onebackground=1;
%else
%    onebackground=0;
% %end
% if ~exist('correct','var') || isempty(correct)  %%MHD
%     correct=[];
% end


% -=- Decide whether avi, stack, or images; set up -=---------------------
[filepath,junk,ext]=fileparts(inputnames);
names=dir(inputnames);
if strcmpi(ext,'.avi')
    movtype='avi';
    v = VideoReader(fullfile(filepath,names.name));
    color_depth=v.BitsPerPixel;
    
    if strcmp(class(read(v,1)),'uint8')
        color_depth=2^8;
    end
    
    ht=v.Height;
    wd=v.Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) v.NumFrames]);
elseif numel(names)==1 && ( strcmpi(ext,'.tif') || ...
        strcmpi(ext,'.tiff') || strcmpi(ext,'.gif') ) % single file, looks like an image stack
    movtype='stack';
    movinfo=imfinfo(fullfile(filepath,names.name));
    color_depth=2^(movinfo(1).BitDepth);
    ht=movinfo(1).Height;
    wd=movinfo(1).Width;
    tmin=max([framerange(1) 1]);
    tmax=min([framerange(2) length(movinfo)]);
elseif numel(names)<6 && ( strcmpi(ext,'.tif') || ...
        strcmpi(ext,'.tiff') || strcmpi(ext,'.gif') ) %MHD 5/21/2018 multiple tiff stacks %single file, looks like an image stack
    movtype='stacks';
    movinfo=imfinfo(fullfile(filepath,names(1).name));
    color_depth=2^(movinfo(1).BitDepth);
    ht=movinfo(1).Height;
    wd=movinfo(1).Width;
    tmin=framerange(1);
    tmax=framerange(2);
    Stacks_lengths=numel(movinfo);
    for jj=2:length(names)
        Stacks_lengths(jj)=length(imfinfo(fullfile(filepath,names(jj).name)));
    end
    Nf=tmax-tmin+1; % frame count
    
    whichstack=ones(Nf,1);
    for jj=2:length(Stacks_lengths)
        whichstack(sum(Stacks_lengths(1:(jj-1)))+1:sum(Stacks_lengths(1:jj)))=jj;
    end
    whichframe=zeros(Nf,1);
    whichframe(1:Stacks_lengths(1))=1:Stacks_lengths(1);
    for jj=2:length(Stacks_lengths)
        whichframe(sum(Stacks_lengths(1:(jj-1)))+1:sum(Stacks_lengths(1:jj)))=...
            1:Stacks_lengths(jj);
    end
    
elseif strcmp(ext,'.seq')
    movtype='sequence';
    
    fid = fopen(inputnames,'r','b'); % open the sequence
    endianType = 'ieee-le'; % use little endian machine format ordering for reading bytes
    
    fseek(fid,548,'bof');  % jump to position 548 from beginning
    imageInfo = fread(fid,24,'uint32',0,endianType);
    headerInfo.ImageWidth = imageInfo(1); % unused
    headerInfo.ImageHeight = imageInfo(2); % unused
    headerInfo.ImageBitDepth = imageInfo(3); % unused
    headerInfo.ImageBitDepthReal = imageInfo(4);
    headerInfo.ImageSizeBytes = imageInfo(5);
    vals = [0,100,101,200:100:600,610,620,700,800,900];
    fmts = {'Unknown','Monochrome','Raw Bayer','BGR','Planar','RGB',...
        'BGRx', 'YUV422', 'YUV422_20', 'YUV422_PPACKED', 'UVY422', 'UVY411', 'UVY444'};
    headerInfo.ImageFormat = fmts{vals == imageInfo(6)};
    fseek(fid,572,'bof'); % headerInfo: get amount of allocated frames
    headerInfo.AllocatedFrames = fread(fid,1,'ushort',endianType);
    fseek(fid,620,'bof'); % headerInfo: get sequence compression
    headerInfo.Compression = fread(fid,1,'uint8',endianType);
    numPixels = headerInfo.ImageSizeBytes;
    color_depth=2^(headerInfo.ImageBitDepthReal);
    ht=headerInfo.ImageHeight;
    wd=headerInfo.ImageWidth;
    switch headerInfo.ImageBitDepthReal % set bit depth for reading
    case 8
        bitDepth = 'uint8';
    case {10,12,14,16}
        bitDepth = 'uint16'; % not tested
    end
    
    numPixels=ht*wd*2;
    
    tt=[1:headerInfo.AllocatedFrames];
    tmin=max([framerange(1) min(tt)]);
    tmax=min([framerange(2) max(tt)]);
    
    fclose(fid);
else
    movtype='images';
    movinfo=imfinfo(fullfile(filepath,names(1).name));
    color_depth=2^(movinfo.BitDepth);
    ht=movinfo.Height;
    wd=movinfo.Width;
    for ii=1:numel(names);
        if strcmp(fullfile(filepath,names(ii).name),bground_name)
            names(ii)=[]; % don't try to track the background file
            break
        end
    end
    tt=NaN(size(names));
    for ii=1:numel(names)
        [junk,myname]=fileparts(names(ii).name);
        tt(ii)=str2double(myname);
    end
    if any(isnan(tt)) % at least one filename was not a number, so just make ordinals
        tt=1:numel(tt);
    end
    tmin=max([framerange(1) min(tt)]);
    tmax=min([framerange(2) max(tt)]);
    names = names( (tt>=tmin) & (tt<=tmax) );
end
Nf=tmax-tmin+1; % frame count

% -=- Pre-compute logarithms for locating particle centers -=-------------
if arealim==1
    logs = 1:color_depth;
    logs = [log(0.0001) log(logs)];
else
    logs=[];
end

% -=- Read background image -=--------------------------------------------
if onebackground %%MHD
    if exist(bground_name,'file')==2
        
        if ~isempty(correct)
            background =double(undistortImage(imread(bground_name),correct));
        else
            background = double(imread(bground_name));
        end
        
    else
        warning('MATLAB:ParticleFinder:noBackgroundFile', ...
            ['Cannot find ' bground_name '. Using blank image instead.'])
        if invert==1
            background = color_depth*ones(ht,wd);
        else
            background = zeros(ht,wd);
        end
    end % if exist(bground_name,'file')==2
    if ndims(background)==3 % check for RGB instead of grayscale
        background=round(mean(background,3));
    end
end
% -=- Set up output file and variables -=---------------------------------
% if writefile
%     if exist(outputname,'file')
%         yn=input(['File ' outputname ' exists. Overwrite (y/n)? '],'s');
%         if ~strcmpi(yn(1),'y')
%             disp('Not overwritten.')
%             writefile=false;
%         else
%             disp(['Replacing file ' outputname '.'])
%             fid=fopen(outputname,'w');
%             fwrite(fid,Nf,'int32'); % header is integer frame count
%         end
%         pause(1)
%     else
%         fid=fopen(outputname,'w');
%         fwrite(fid,Nf,'int32'); % header is integer frame count
%     end
% end
%begins=ones(Nf+1,1);
%mark=[];

%MHD do in order of each background image.

% inds=[];
% for i=1:max(index)
%     inds=[inds, find(index==i)];
% end

f=0; %frame

for ii=1:Nf
    f=f+1; %frame
    tt=tmin+ii-1; % current time
    switch movtype
        case('avi')
            
            im = double(read(v,tt));
            
            NAME = names.name;
        case('stack')
            if ~isempty(correct)
                im = double(undistortImage(imread(fullfile(filepath,names.name),tt),correct));
                
            else
                im = double(imread(fullfile(filepath,names.name),tt));
            end
            
            imoriginal = imread(fullfile(filepath,names.name),tt);
            NAME= names.name;
        case('stacks') %MHD
            k=whichstack(ii);
            m=whichframe(ii);
            
            %im = double(imread(fullfile(filepath,names(k).name),m));
            
            if ~isempty(correct)
                im =double(undistortImage(imread(fullfile(filepath,names(k).name),m),correct));
            else
                im = double(imread(fullfile(filepath,names(k).name),m));
            end
            
            imoriginal = (imread(fullfile(filepath,names(k).name),m));
            
            NAME=[names(k).name '-' num2str(m)];
            
        case('sequence')
            
            [imoriginal,~]=loadseqframe(inputnames,ii);
            im=double(imoriginal);
%                   first=1; %???? mhd  
%         % set readStart of next frame
%         readStart = ((first-1+ii-1)*numPixels)+((first-1+ii)*8192);
%         
%         % read uncompressed image
%         fseek(fid,readStart,'bof');
%         UYVYcolVec = fread(fid, numPixels/2, bitDepth, endianType); %mimi changed numpixels/2
%        
%         
%                bd=headerInfo.ImageBitDepthReal;
%                
%                keyboard
% 
%         if length(UYVYcolVec)==ht*wd
%     
%                 if bd==8
%                     im = uint8(reshape(UYVYcolVec,headerInfo.ImageWidth,headerInfo.ImageHeight)');
%                 else
%                     im = uint16(reshape(UYVYcolVec,headerInfo.ImageWidth,headerInfo.ImageHeight)');
%                 end
%                 
%                 imoriginal=im;
%                 im=double(im);
%             elseif length(UYVYcolVec)==ht*wd+12
%                 
%                 if bd==8
%                     im = uint8(reshape(UYVYcolVec(1:ht*wd),headerInfo.ImageWidth,headerInfo.ImageHeight)');
%                 else
%                     
%                     im = uint16(reshape(UYVYcolVec(1:ht*wd),headerInfo.ImageWidth,headerInfo.ImageHeight)');
%                 end
%                   imoriginal=im;
%                 im=double(im);
%         
%         
%         else
%             im=zeros(ht,wd);
% 
%         end
        %SEQ = uint8(UYVYimg); %mimi;
% %         SEQ=UYVYimg;
        
        % read timestamp
        %readStart = readStart + numPixels;
        %fseek(fid,readStart,'bof');
       % time = readTimestamp(fid, endianType);
       % I{i,2} = time;
%        im = doubleSEQ;

        NAME='seq';
            
        case('images')
            %im = double(imread(fullfile(filepath,names(ii).name)));
            
            %MHD
            if ~isempty(correct)
                im =double(undistortImage(imread(fullfile(filepath,names(ii).name)),correct));
            else
                im = double(imread(fullfile(filepath,names(ii).name)));
            end
            
            imoriginal = (imread(fullfile(filepath,names(ii).name)));
            
            NAME=names(ii).name;
            
    end % switch movtype
    if ndims(im)==3
        if size(im,3)>3
            error('Sorry, only grayscale and RGB images are supported.')
        end
        im=round(mean(im,3)); % convert to grayscale if necessary
    end
    if invert==1 && logical(onebackground) % dark particles on light background
        im = background - im;
    elseif invert==0 && logical(onebackground)% light particles on dark background
        im = im - background;
    elseif onebackground==0  %phase averaged background- MHD
        
        ind=index(ii);  %phase
        if ii==1 || (ind~=index(lastind)) %load new background image
            
            bgname=[bground_name num2str(ind) '.tif'];
            
            %%MHD
            if exist(bgname,'file')==2
                %MHD
                if ~isempty(correct)
                    background =double(undistortImage(imread(bgname),correct));
                else
                    background = double(imread(bgname));
                end
                
                
            else
                warning('MATLAB:ParticleFinder:noBackgroundFile', ...
                    ['Cannot find ' bground_name '. Using blank image instead.'])
                if invert==1
                    background = color_depth*ones(ht,wd);
                else
                    background = zeros(ht,wd);
                end
            end % if exist(bground_name,'file')==2
            if ndims(background)==3 % check for RGB instead of grayscale
                background=round(mean(background,3));
            end
        end
        
        if invert==1
            im = background - im;
        elseif invert==0
            im = im - background;
        end
    else
        im = abs( im - background ); % seek any contrast, light or dark
    end
    im(im<0) = 0;
    
 
    
    
    %M&K
    if arealim>1
        [pos,ang1] = FindRegions(im,threshold,arealim);
    else
        pos = FindParticles(im,threshold,logs);
        ang1 = [];
    end
    
    
    %ORIGINAL
    %     if arealim==1
    %         pos = FindParticles(im,threshold,logs);
    %     else
    %         [pos,ang1] = FindRegions(im,threshold,arealim); % could also keep angle...
    %     end
    
    
    N=size(pos,1);
    if ii==1 % if first frame, pre-allocate arrays for speed
        x=NaN(N*Nf,1);
        y=NaN(N*Nf,1);
        t=NaN(N*Nf,1);
        if arealim~=1
            ang=NaN(N*Nf,9); %MHD 5/16/2018 added perimeter and weightedcentroid outputs
        else
            ang=[]; % no angles
        end
        memloc=1;
    end
    if N>0
        x(memloc:memloc+N-1)=pos(:,1);
        y(memloc:memloc+N-1)=pos(:,2);
        t(memloc:memloc+N-1)=tt;
        if arealim~=1  %%MHD commented out
            ang(memloc:memloc+N-1,:)=ang1;
        end
        memloc=memloc+N;
    end
    % begins(ii+1)=begins(ii)+N;
    % display progress
    if rem(f,25)==0 %display every 25 frames.
        %         switch movtype
        %             case('avi')
        %                 disp(['Found ' num2str(N,'%.0f') ' particles in ' ...
        %                     names.name ' frame ' num2str(tt) ' ('  ...
        %                     num2str(ii) ' of ' num2str(Nf) ').'])
        %             case('stack')
        %                 disp(['Found ' num2str(N,'%.0f') ' particles in ' ...
        %                     names.name ' frame ' num2str(tt) ' (' ...
        %                     num2str(ii) ' of ' num2str(Nf) ').'])
        %             case('images')
        %                 disp(['Found ' num2str(N,'%.0f') ' particles in ' ...
        %                     names(ii).name ' (' num2str(ii) ' of ' num2str(Nf) ').'])
        %         end
        
        
        %MHD
        disp(['Found ' num2str(N,'%.0f') ' particles in ' NAME ' frame (' num2str(f) ' of ' num2str(Nf) ').'])
        
    end% switch movtype
    %     if writefile
    %         fwrite(fid,N,'int32'); % start the frame w/ particle count
    %         fwrite(fid,pos','float32'); % then write particle locations x1 y1 x2 y2  ...
    %         % Next: write angle as well!
    %     end
    lastind=ii;
end % for ii=1:tmax
%fclose(fid);

x(memloc:end)=[];
y(memloc:end)=[];
t(memloc:end)=[];
if arealim~=1
    ang(memloc:end,:)=[];
end

% MHD put back in chronological order
[~,I]=sort(t);
x=x(I);
y=y(I);
t=t(I);
% LKC if statement so no angle for tracers
%(I THINK??)
if arealim ~= 1
    ang=ang(I,:);
else
    ang = [];
end


% if writefile
%     fclose(fid);
% end
% 
% % -=- Plot if requested -=------------------------------------------------
% if noisy
%     if isnumeric(noisy) && noisy>1
%         disp(['Plotting and saving frames. ' ...
%             'Please do not cover the figure window!'])
%         if exist(savedirname,'file')~=7
%             mkdir(savedirname)
%         end
%     else
%         disp('Plotting...')
%     end
%     figure;
%     axes('nextplot','add','dataaspectratio',[1 1 1],'ydir','reverse', ...
%         'xlim',0.5+[0 size(im,2)],'ylim',0.5+[0 size(im,1)]);
%     switch movtype
%         case('avi')
%             xlabel([names.name ', threshold = ' num2str(threshold) ...
%                 ', ' bground_name ', arealim = ' num2str(arealim), ...
%                 ', invert = ' num2str(invert)],'interpreter','none');
%             [im,mark] = read_uncompressed_avi( ...
%                 fullfile(filepath,names.name),1);
%             im = im.cdata;
%         case('stack')
%             xlabel([names.name ', threshold = ' num2str(threshold) ...
%                 ', ' bground_name ', arealim = [' num2str(arealim), ...
%                 '], invert = ' num2str(invert)],'interpreter','none');
%             im = imread(fullfile(filepath,names.name),1);
%         case('images')
%             xlabel([inputnames ', threshold = ' num2str(threshold) ...
%                 ', ' bground_name ', arealim = [' num2str(arealim) ']'], ...
%                 'interpreter','none');
%             im = imread(fullfile(filepath,names(1).name));
%     end % switch movtype
%     hi=imagesc(im);
%     colormap(gray); % has no effect if image is color b/c RGB values override
%     hl=plot(NaN,NaN,'r.');
%     mark=[];
%     for ii=1:Nf
%         ind=begins(ii):begins(ii+1)-1;
%         if arealim==1 % just dots, no angles
%             set(hl,'xdata',x(ind),'ydata',y(ind));
%         else % plot bars to show angles
%             delete(hl);
%             hl=NaN(size(ind));
%             for jj=1:numel(hl)
%                 hl(jj)=plot(x(ind(jj))+[-1 1]*barsize*cos(ang(ind(jj))), ...
%                     y(ind(jj))+[-1 1]*barsize*sin(ang(ind(jj))),'r-');
%             end
%         end % if arealim==1
%         switch movtype
%             case('avi')
%                 [im,mark] = read_uncompressed_avi( ...
%                     fullfile(filepath,names.name),tmin+ii-1,mark);
%                 set(hi,'cdata',im.cdata);
%                 set(gcf,'name',[num2str(numel(ind)) ' particles in ' ...
%                     names.name ' frame ' num2str(tmin+ii-1) ' (' ...
%                     num2str(ii) ' of ' num2str(Nf) ')']);
%             case('stack')
%                 set(hi,'cdata',imread(fullfile(filepath,names.name), ...
%                     tmin+ii-1));
%                 set(gcf,'name',[num2str(numel(ind)) ' particles in ' ...
%                     names.name ' frame ' num2str(tmin+ii-1) ' (' ...
%                     num2str(ii) ' of ' num2str(Nf) ')']);
%             case('images')
%                 set(hi,'cdata',imread(fullfile(filepath, ...
%                     names(ii).name)));
%                 set(gcf,'name',[num2str(numel(ind)) ' particles in ' ...
%                     names(ii).name ' (' num2str(ii) ' of ' ...
%                     num2str(Nf) ')']);
%         end % switch movtype
%         drawnow
%         pause(pausetime);
%         if isnumeric(noisy) && noisy>1
%             snap=getframe(gca);
%             imwrite(snap.cdata, ...
%                 fullfile(filepath,savedirname,names(ii).name));
%         end
%     end
% end
disp('Done.')

end % function ParticleFinder

% -=- function FindParticles -=-------------------------------------------
function pos = FindParticles(im, threshold, logs)
% Given an image "im", FindParticles finds small particles that are
% brighter than their four nearest neighbors and also brighter than
% "threshold". Particles are located to sub-pixel accuracy by applying a
% Gaussian fit in each spatial direction. The input "logs" depends on
% the color depth and is re-used for speed. Particle locations are
% returned in the two-column array "pos" (with x-coordinates in the first
% column and y-coordinates in the second).

s = size(im);

% identify the local maxima that are above threshold
maxes = find(im >= threshold & ...
    im > circshift(im,[0 1]) & ...
    im > circshift(im,[0 -1]) & ...
    im > circshift(im,[1 0]) & ...
    im > circshift(im,[-1 0]));

% now turn these into subscripts
[x,y] = ind2sub(s, maxes);

% throw out unreliable maxes in the outer ring
good = find(x~=1 & y~=1 & x~=s(1) & y~=s(2));
x = x(good);
y = y(good);

% find the horizontal positions

% look up the logarithms of the relevant image intensities
z1 = logs(im(sub2ind(s,x-1,y)) + 1)';
z2 = logs(im(sub2ind(s,x,y)) + 1)';
z3 = logs(im(sub2ind(s,x+1,y)) + 1)';

% compute the centers
xcenters = -0.5 * (z1.*(-2*x-1) + z2.*(4*x) + z3.*(-2*x+1)) ./ ...
    (z1 + z3 - 2*z2);

% do the same for the vertical position
z1 = logs(im(sub2ind(s,x,y-1)) + 1)';
z3 = logs(im(sub2ind(s,x,y+1)) + 1)';
ycenters = -0.5 * (z1.*(-2*y-1) + z2.*(4*y) + z3.*(-2*y+1)) ./ ...
    (z1 + z3 - 2*z2);

% make sure we have no bad points
good = find(isfinite(xcenters) & isfinite(ycenters));

% fix up the funny coordinate system used by matlab
pos = [ycenters(good), xcenters(good)];

end % function FindParticles


%% michael & kasper's code
function [pos,ang] = FindRegions(im,threshold,arealim)
% Given an image "im", FindRegions finds regions that are brighter than
% "thresold" and have area larger than "arealim". Region centroids are
% returned in the two-column array "pos" (with x-coordinates in the first
% column and y-coordinates in the second). Region orientations are
% returned in radians, in the vector "ang".
debug=0;
%blobyes; %get blob props



% imagesc(im)
% pause


if numel(arealim)==1
    arealim=[arealim inf]; % assume single size is a minimum  %CAN inPUT area min & max
end
s = size(im);
warnstate=warning('off','MATLAB:divideByZero'); % squelch regionprops divide-by-zero warnings
inm = im>threshold;
inm = bwareaopen(inm,min(arealim));
%
%     if areaweighting234
%
%     props=regionprops(inm,im,'WeightedCentroid','Area', ...
%         'PixelValues','PixelList','Orientation','MajorAxisLength','MinorAxisLength');
%
%
%      pos=vertcat(props.WeightedCentroid);
%
%  else


%%%%%% MHD 08/03/2023 --- adding option for weighted centroid tracking
props=regionprops(inm,im,'Centroid','Area', ...
    'MeanIntensity','MaxIntensity','Orientation','MajorAxisLength','MinorAxisLength','WeightedCentroid','Perimeter');

weightedcentroid = 1;
if weightedcentroid

    % this is new option
pos=vertcat(props.WeightedCentroid);

else
    % this is original tracking code: uses normal centroid
pos=vertcat(props.Centroid);

end


%MHD
ori= cat(1,props.Orientation);
a  = cat(1,props.MajorAxisLength);
b  = cat(1,props.MinorAxisLength);
are=cat(1,props.Area);
%MHD added perimeter and wc outputs 5/16/2018
peri=cat(1,props.Perimeter);
wc=vertcat(props.WeightedCentroid);
mi=cat(1,props.MeanIntensity);
maxi=cat(1,props.MaxIntensity);
if numel(pos)>0  %are there any found regions?
    good = pos(:,1)~=1 & pos(:,2)~=1 & pos(:,1)~=s(2) & ...
        pos(:,2)~=s(1) & vertcat(props.Area)>arealim(1) & ...
        vertcat(props.Area)<arealim(2); % remove regions on edge, too small, or too big
    pos     = pos(good,:);
    %props   = props(good);
    ori     = ori(good);
    a       = a(good);
    b       = b(good);
    are     = are(good);
    peri    = peri(good);
    wc      = wc(good,:);
    mi      = mi(good);
    maxi    = maxi(good);
    ang     = [ori a b are peri mi maxi wc];
else
    ang = [];
end



% show for debuggin
if debug
    % test = 1
    %   fh= figure(1);clf(1);imagesc(im); hold on;
    fh= figure(1);clf(1);imagesc(imoriginal); hold on; %or inm
    set(fh,'Position',[500 100 700 700]);
    if size(pos,1)>0
        plot(pos(:,1),pos(:,2),'ro');   %particle positions
    end
    %  set(gca,'Clim',[0 50]);
    magg = 20;
    hq = quiver(pos(:,1),pos(:,2),magg*cosd(ori(:)),-magg*sind(ori(:)),0);
    set(hq,'color','m');
    keyboard
end
%
%
warning(warnstate) % put things back the way we found them
end % function FindRegions




% -=- function FindRegions -=---------------------------------------------
% % function [pos,ang] = FindRegions(im,threshold,arealim)
% % % Given an image "im", FindRegions finds regions that are brighter than
% % % "thresold" and have area larger than "arealim". Region centroids are
% % % returned in the two-column array "pos" (with x-coordinates in the first
% % % column and y-coordinates in the second). Region orientations are
% % % returned in radians, in the vector "ang".
% %
% %     if numel(arealim)==1
% %         arealim=[arealim inf]; % assume single size is a minimum
% %     end
% %
% %     s = size(im);
% %     warnstate=warning('off','MATLAB:divideByZero'); % squelch regionprops divide-by-zero warnings
% %     props=regionprops(im>threshold,im,'WeightedCentroid','Area', ...
% %         'PixelValues','PixelList');
% %     pos=vertcat(props.WeightedCentroid);
% %     if numel(pos)>0
% %         good = pos(:,1)~=1 & pos(:,2)~=1 & pos(:,1)~=s(2) & ...
% %             pos(:,2)~=s(1) & vertcat(props.Area)>arealim(1) & ...
% %             vertcat(props.Area)<arealim(2); % remove regions on edge, too small, or too big
% %         pos=pos(good,:);
% %         props=props(good);
% %     end
% %     ang = NaN(size(props));
% %     for ii = 1:numel(ang)
% %         xy = props(ii).PixelList;
% %         vals = props(ii).PixelValues;
% %         uxx = sum( (xy(:,1)-pos(ii,1)).^2 .* vals ,1) / sum(vals) + 1/12;
% %         uyy = sum( (xy(:,2)-pos(ii,2)).^2 .* vals ,1) / sum(vals) + 1/12;
% %         uxy = sum( (xy(:,1)-pos(ii,1)).*(xy(:,2)-pos(ii,2)) .* vals ,1) ...
% %             / sum(vals);
% %         if (uyy > uxx)
% %             num = uyy - uxx + sqrt((uyy - uxx)^2 + 4*uxy^2);
% %             den = 2*uxy;
% %         else
% %             num = 2*uxy;
% %             den = uxx - uyy + sqrt((uxx - uyy)^2 + 4*uxy^2);
% %         end
% %         if (num == 0) && (den == 0)
% %             ang(ii) = 0;
% %         else
% %             ang(ii) = atan(num/den);
% %         end
% %     end % for ii = 1:numel(ang)
% %     warning(warnstate) % put things back the way we found them
% % end % function FindRegions

