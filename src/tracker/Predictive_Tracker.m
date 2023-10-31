function [vtracks,ntracks,meanlength,rmslength] = Predictive_Tracker(centroids,max_disp)
%% things to work on
% separating tracking and finding and debugging a little better
% more analysis code
% more flexible inputs, could even include json inputs
% more options for using e.g., centroid, weighted centroid, PIV, frame
% straddling

%% info
% Predictive tracker in 2D given identified particle locations

% Inputs:

% centroids should be a structure where
% centroids.x is the x location vector
% centroids.y is the y location vector
% centroids.t is the frame number vector (should be integers)

% centroids can have other info that will be passed to the output tracks e.g.,
% centroids.area

% everything should be column vectors/matrices with the same column length.
% the extra fields can have different row lengths

%% other inputs to parse

fitwidth = 3; % 3 nick default
filterwidth = 1;
yesvels = 1;

%% Parse input particle info
try
    x=centroids.x;
    y=centroids.y;
    t=centroids.t;

    % Make sure at least x,y,t are column vectors
    if ~iscolumn(x);  x=x'; end

    if ~iscolumn(y);  y=y'; end

    if ~iscolumn(t);  t=t'; end

catch
    error('input centroids structure is not properly formatted. needs to include ''x'', ''y'', and ''t''')
end

% Make sure they are all the same size
N_found = length(x);

if N_found ~= length(y) || N_found~= length(t)
    error('input ''x'', ''y'', and ''t'' are not all the same length')
elseif sum(rem(t,1))>0
    error('input time steps are not integer frame numbers')
end

% sort anything else included in centroids
centroid_fields_all = fieldnames(centroids);
centroid_fields_extra = {};



if length(centroid_fields_all) > 3
    for i=1:length(centroid_fields_all)

        field_i = centroid_fields_all{i};

        if ~sum(strcmp(field_i,{'x','y','t'}))
            field_data = getfield(centroids,field_i);

            [m,n] = size(field_data);

            if m == N_found
                centroid_fields_extra = [centroid_fields_extra, field_i];

            elseif n == N_found %need to transpose the data

                warning(['transposed input data in ' field_i])

                centroid_fields_extra = [centroid_fields_extra, field_i];
                centroids = setfield(centroids,field_i, field_data');
            else
                warning(['ignoring input data in ' field_i ' because it is not the right size'])
            end
        end
    end
end
N_extra = length(centroid_fields_extra);



% Make sure everything is sorted in increasing time
if ~issorted(t)
    [t,ind] = sort(t);
    x = x(ind);
    y = y(ind);
    for f = 1:N_extra
        field_f = centroid_fields_extra{f};
        field_data = getfield(centroids,field_f);
        field_data = field_data(ind,:);
        centroids = setfield(centroids,field_f, field_data);
    end
end

    %% Sort out which particles were found in which frames

    % Start frame counter at 1 (helps streamline code)
    start_time = min(t) -1;
    t = t - start_time;

    tt2=1:max(t);

    [~,endsind]=unique(t,'last');
    [tt,beginsind]=unique(t,'first');

    ends=zeros(length(tt2),1);
    begins=ends;

    ends(tt)=endsind;
    begins(tt)=beginsind;
    tt=tt2;

    Nf=numel(tt);
    if Nf < (2*fitwidth+1)
        error(['Sorry, found too few files named ' inputnames '.'])
    end


    %% Set up array for tracks using particles in the first frame

    extra_inputs_init = cell(1, 2*N_extra);
    extra_inputs_init(1,1:2:end) = centroid_fields_extra;

    % find particles in first frame
    ind = begins(1):ends(1);
    npart_frame = numel(ind);  % originally nparticles in nick code

    % start new tracks with particles from first frame
    tracks = repmat(struct('len',[],'X',[],'Y',[],'T',[], extra_inputs_init{:}), ...
        npart_frame,1);

    for ii = 1:npart_frame
        extra_inputs = struct_parser(centroid_fields_extra, centroids, ind(ii));
        tracks(ii) = struct('len',1,'X',x(ind(ii)),'Y',y(ind(ii)),'T',1, ...
            extra_inputs{:});
    end


    %% -=- Keep track of which tracks are active -=----------------------------
    active = 1:npart_frame;
    n_active = numel(active);
    disp(['Processed frame 1 of ' num2str(Nf) '.'])
    disp(['    Number of particles found: ' num2str(npart_frame,'%.0f')])
    disp(['    Number of active tracks: ' num2str(n_active,'%.0f')])
    disp(['    Total number of tracks: ' num2str(numel(tracks),'%.0f')])

    %% -=- Loop over frames, starting in second frame -=-----------
    for t = 2:Nf

        ind=begins(t):ends(t);
        time=tt(t);

        if begins(t)==1 & t~=1
            nfr1=0;
            ind=[];
        else
            nfr1 = numel(ind);
        end

        if nfr1==0
            warning('MATLAB:PredictiveTracker:noParticles', ...
                ['Found no particles in frame ' num2str(t) '.']);
        end

        fr1=[x(ind) y(ind)];


        centroid_extras = struct();
        for f = 1:N_extra
            field_data = getfield(centroids,centroid_fields_extra{f});
            centroid_extras = setfield(centroid_extras,centroid_fields_extra{f},field_data(ind,:));
        end


        %% -=- Match the tracks with kinematic predictions -=----------------------

        % for convenience, we'll grab the relevant positions from the tracks
        now = zeros(n_active,2);
        prior = zeros(n_active,2);
        for ii = 1:n_active
            tr = tracks(active(ii));
            now(ii,1) = tr.X(end);
            now(ii,2) = tr.Y(end);
            if tr.len > 1
                prior(ii,1) = tr.X(end-1);
                prior(ii,2) = tr.Y(end-1);
            else
                prior(ii,:) = now(ii,:);
            end
        end

        % estimate a velocity for each particle in fr0
        velocity =now - prior;  %Use this for continuous

        estimate = now + velocity;

        % define cost and link arrays
        costs = zeros(n_active,1);
        links = zeros(n_active,1);

        % keyboard
        if nfr1>0
            % loop over active tracks
            for ii = 1:n_active
                % now, compare this estimated positions with particles in fr1
                dist_fr1 = (estimate(ii,1)-fr1(:,1)).^2 + (estimate(ii,2)-fr1(:,2)).^2;
                % save its cost and best match
                costs(ii) = min(dist_fr1);
                if costs(ii) > max_disp^2
                    continue;
                end
                bestmatch = find(dist_fr1 == costs(ii));
                % if there is more than one best match, we are confused; stop
                if numel(bestmatch) ~= 1
                    continue;
                end
                % has another track already matched to this particle?
                ind = links == bestmatch;
                if sum(ind) ~= 0
                    if costs(ind) > costs(ii)
                        % this match is better
                        links(ind) = 0;
                    else
                        continue;
                    end
                end
                links(ii) = bestmatch;
            end

            % now attach the matched particles to their tracks
            matched = zeros(nfr1,1);
            for ii = 1:n_active
                if links(ii) ~= 0
                    % this track found a match
                    tracks(active(ii)).X(end+1) = fr1(links(ii),1);
                    tracks(active(ii)).Y(end+1) = fr1(links(ii),2);
                    tracks(active(ii)).len = tracks(active(ii)).len + 1;
                    tracks(active(ii)).T(end+1) = time;

                    for f = 1:N_extra
                        field_data = getfield(centroid_extras,centroid_fields_extra{f});
                        tracks(active(ii)) = setfield(tracks(active(ii)),centroid_fields_extra{f},[getfield(tracks(active(ii)),centroid_fields_extra{f}); field_data(links(ii),:)]);
                    end

                    matched(links(ii)) = 1;
                end
            end
            active = active(links~=0);

            % and start new tracks with the particles in fr1 that found no match
            unmatched = find(matched == 0);

            newtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[], ...
                extra_inputs_init{:}),numel(unmatched),1);
            for ii = 1:numel(unmatched)
                extra_inputs = struct_parser(centroid_fields_extra, centroids, unmatched(ii));
                newtracks(ii) = struct('len',1,'X',fr1(unmatched(ii),1),...
                    'Y',fr1(unmatched(ii),2),'T',time,extra_inputs{:});
            end
            %end
        else % if nfr1>0
            active=[];
            newtracks=[];
            unmatched=[];
        end
        active = [active (numel(tracks)+1):(numel(tracks)+numel(newtracks))];
        tracks = [tracks ; newtracks];
        n_active = numel(active);

        disp(['Processed frame ' num2str(t) ' of ' num2str(Nf) '.'])
        disp(['    Number of particles found: ' num2str(nfr1,'%.0f')])
        disp(['    Number of active tracks: ' num2str(n_active,'%.0f')])
        disp(['    Number of new tracks started here: ' ...
            num2str(numel(unmatched),'%.0f')])
        disp(['    Number of tracks that found no match: ' ...
            num2str(sum(links==0),'%.0f')])
        disp(['    Total number of tracks: ' num2str(numel(tracks),'%.0f')])



    end

    if ~yesvels
        vtracks = tracks([tracks.len] > 1); %remove unmatched
        ntracks=numel(vtracks);
        meanlength = mean([vtracks.len]);
        rmslength = sqrt(mean([vtracks.len].^2));

        %%% look for errors

        %     for kk=1:length(vtracks)
        %         T=vtracks(kk).T;
        %        if sum(diff(T)~=1)>0
        %            keyboard
        %        end
        %
        %        XX=vtracks(kk).X;
        %        if XX(2)==XX(1)
        %            keyboard
        %        end
        %
        %     end


    else

        % -=- Prune tracks that are too short -=----------------------------------
        disp('Pruning...');
        tracks = tracks([tracks.len] >= (2*fitwidth+1));
        ntracks = numel(tracks);
        meanlength = mean([tracks.len]);
        rmslength = sqrt(mean([tracks.len].^2));

        % -=- Compute velocities -=-----------------------------------------------
        disp('Differentiating...');

        % define the convolution kernel
        Av = 1.0/(0.5*filterwidth^2 * ...
            (sqrt(pi)*filterwidth*erf(fitwidth/filterwidth) - ...
            2*fitwidth*exp(-fitwidth^2/filterwidth^2)));
        vkernel = -fitwidth:fitwidth;
        vkernel = Av.*vkernel.*exp(-vkernel.^2./filterwidth^2);

        % loop over tracks
        %if minarea==1
        vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[]),ntracks,1);
        %else
        vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[], ...
            'V',[],extra_inputs_init{:}),ntracks,1);
        %end
        for ii = 1:ntracks
            u = -conv(tracks(ii).X,vkernel,'valid');
            v = -conv(tracks(ii).Y,vkernel,'valid');


            ind = (fitwidth+1):(tracks(ii).len-fitwidth);

            extra_inputs = struct_parser(centroid_fields_extra, tracks(ii),ind);

            vtracks(ii) = struct('len',tracks(ii).len - 2*fitwidth, ...
                'X',tracks(ii).X(ind), ...
                'Y',tracks(ii).Y(ind), ...
                'T',tracks(ii).T(ind), ...
                'U',u, ...
                'V',v, extra_inputs{:});
        end
    end


    %vtracks = tracks;
    %vtracks(:).T = vtracks(:).T+start_time;
    % add start_time back to particle time at the end
end


    function [extra_inputs] = struct_parser(fnames, data, ind)
        % make command to append extra data to output tracks
        % input is centroid_fields_extra and index in string
        % data is centroids

        extra_inputs = cell(1, 2*length(fnames));
        extra_inputs(1,1:2:end) = fnames;

        for i = 1:length(fnames)
            field_data = getfield(data, fnames{i});
            extra_inputs(1,i*2) = {field_data(ind,:)};
        end


    end