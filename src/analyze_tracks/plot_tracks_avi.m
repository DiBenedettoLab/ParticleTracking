function [] = plot_tracks_avi(inputname,vtracks,framerange)
%framerange is vector e.g. [5 10] which denotes first and last frame to
%plot


vid = VideoReader(inputname);
if isempty(framerange)
   t = 1:vid.Numframes;
else
   t = framerange(1):framerange(end);
end

colorlist=get(gca,'colororder');
Ncolors=size(colorlist,1);

if ~isempty(vtracks)
    ntracks=length(vtracks);
    framerange=NaN(ntracks,2);
    for jj=1:ntracks
        framerange(jj,1)=vtracks(jj).T(1);
        framerange(jj,2)=vtracks(jj).T(end);
    end
end

first_track = min(framerange(:));
last_track  = max(framerange(:));


for i=t
    
    if i<first_track
        continue
    elseif i > last_track
        break
    end
    
    im = read(vid,i);
        
       
        imagesc(im)
        colormap(gray)
        
        
   
        hold on
        ind = find( (framerange(:,1)<=i) & (framerange(:,2)>=i) );
        for jj=1:numel(ind)
            col=colorlist(mod(ind(jj)-1,Ncolors)+1,:);
            indt=1:(i-framerange(ind(jj),1)+1); % plot from beginning of track to current frame
            plot(vtracks(ind(jj)).X(indt),(vtracks(ind(jj)).Y(indt)), ...
                '-','color',col);
            
        end
    
    title(['frame ' num2str(i)],'Interpreter','none')
    axis image
  
    drawnow
    pause(.05)
    

    
end


end
