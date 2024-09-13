function [u_corr,u_guess,corr_peaks] = Particle_Velocity_Cross_Corr(x,y,im1,im2,windowsize,medfilt,yesplot)
%cross correlation of particle image to get velocities
%median filter the image first

%inputs, [x,y], vectors of larvae trajectory centroids
% [t], vector of particle trajectory timestamps
%[seqname], name of the sequence of raw images
%[df], space between cross correlated frames (default is 1)
% [windowsize], pixels


%% parse inputs
%windowsize=45; %pixels ?
%medfilt=[5 5]; %size of median filter, pixels. make 0 to remove this step

medfilt = [medfilt(1) medfilt(1)];
winsize=ceil(windowsize/2);
winsize2=winsize;

%for debugging purposes
%yesplot=1;


% frames
next = 2;
first = 1;

%% Guess the velocity
    u_guess=[x(next)-x(first),y(next)-y(first)];

    %% Extract bounding boxes of particle images for cross-correlating
    yy=round(x(first));
    xx=round(y(first));

    [xpix,ypix]=size(im1);

    indx1=[max(xx-winsize,1):min(xx+winsize,xpix)];
    indy1=[max(yy-winsize,1):min(yy+winsize,ypix)];
    part_image1=im1(indx1,indy1);

    yy=round(x(next));
    xx=round(y(next));

    indx2=[max(xx-winsize2,1):min(xx+winsize2,xpix)];
    indy2=[max(yy-winsize2,1):min(yy+winsize2,ypix)];
    part_image2=im2(indx2,indy2);

    %% Get the particle images ready to cross-correlate
    % Other options?
    % remove saturation
    % remove mean
    % remove median

    % Remove mean
    part_image1=part_image1-mean(part_image1(:));
    part_image2=part_image2-mean(part_image2(:));


    % median filter (to remove noise)

    if medfilt(1)>0
        part_image1=medfilt2(part_image1,medfilt);
        part_image2=medfilt2(part_image2,medfilt);
    end

    %% Cross correlate
    ccr=xcorr2(part_image1,part_image2);

    [peak,ind]=max(ccr(:));

    corr_peaks=peak;
    [r,c]=ind2sub(size(ccr),ind);


    %% subpixel fitting: just linear in each direction for now, 
    % this should be updated to be actually 2d and use more info

    % subpixel fitting: rows
    z1 = log(ccr(r-1,c));
    z2 = log(ccr(r,c));
    z3 = log(ccr(r+1,c));

    % compute the centers
    rcent = -0.5 * (z1.*(-2*r-1) + z2.*(4*r) + z3.*(-2*r+1)) ./ ...
        (z1 + z3 - 2*z2);

    % do the same for the colums
    z1 = log(ccr(r,c-1));
    z3 = log(ccr(r,c+1));
    ccent = -0.5 * (z1.*(-2*c-1) + z2.*(4*c) + z3.*(-2*c+1)) ./ ...
        (z1 + z3 - 2*z2);

    


%% alt subpixel fitting
% https://www.mathworks.com/matlabcentral/fileexchange/26504-sub-sample-peak-fitting-2d
K = (ccr(r-1:r+1,c-1:c+1));
% approximate polynomial parameter
A = (K(2,1)+K(1,1)-2*K(1,2)+K(1,3)-2*K(3,2)-2*K(2,2)+K(2,3)+K(3,1)+K(3,3));
B = (K(3,3)+K(1,1)-K(1,3)-K(3,1));
C = (-K(1,1)+K(1,3)-K(2,1)+K(2,3)-K(3,1)+K(3,3));
%d = (2*K(2,1)-K(1,1)+2*K(1,2)-K(1,3)+2*K(3,2)+5*K(2,2)+2*K(2,3)-K(3,1)-K(3,3));
E = (-2*K(2,1)+K(1,1)+K(1,2)+K(1,3)+K(3,2)-2*K(2,2)-2*K(2,3)+K(3,1)+K(3,3));
F = (-K(1,1)-K(1,2)-K(1,3)+K(3,1)+K(3,2)+K(3,3));
% (ys,xs) is subpixel shift of peak location relative to point (2,2)
rcent = (6*B*C-8*A*F)/(16*E*A-9*B^2) + r;
ccent = (6*B*F-8*E*C)/(16*E*A-9*B^2) + c;



%%
    u_corr=[indy2(1)-indy1(1)-(ccent-size(part_image2,2)),indx2(1)-indx1(1)-(rcent-size(part_image2,1))];

    
    j = 1;
    if yesplot
        subplot(1,2,1)
        hold off
        imagesc(indy1,indx1,part_image1)
        %axis([ min([indy1 indy2]) max([indy1 indy2]) min([indx1 indx2]) max([indx1 indx2])])
        hold on


        scatter(x(j),y(j),'y*')
        quiver(x(j),y(j),u_corr(j,1),u_corr(j,2))

        axis image
        axis([ min([indy1 indy2]) max([indy1 indy2]) min([indx1 indx2]) max([indx1 indx2])])

        subplot(1,2,2)
        hold off
        imagesc(indy2,indx2,part_image2)
        axis image
        axis([ min([indy1 indy2]) max([indy1 indy2]) min([indx1 indx2]) max([indx1 indx2])])
        hold on
        scatter(x(next),y(next),'y*')

        colormap(copper)


        
    end

end




