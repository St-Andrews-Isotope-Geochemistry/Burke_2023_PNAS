%% This code uses sulfur isotope data from ice cores to solve for the fraction of stratospheric sulfate in a sulfate peak
%% Written by Andrea Burke
%% When using cite Burke et al. (2023) "High sensitivity of summer temperatures to stratospheric sulfur
%% loading from volcanoes in the Northern Hemisphere." Proceedings of the National Academy of Sciences (PNAS).


%% Read in data table
numVars = 14; % number of columns of data to read in
varNames = {'Core','Eruption','Type','BotDepth','TopDepth','Age','Volume', 'Conc','d34S', 'd34S_err','d33S', 'd33S_err','D33S', 'D33S_err'} ;
varTypes = {'char', 'char', 'char', 'double','double','double','double','double','double','double','double','double','double','double'};
data_range = 'A4:N122';
opts = spreadsheetImportOptions('NumVariables',numVars,...
    'VariableNames',varNames,...
    'VariableTypes', varTypes,...
    'DataRange', data_range);

imported_data = readtable('Burke_2023_PNAS.xlsx', opts);

%% Choose which event(s) to analyze and set up Monte Carlo simulation parameters

cores = {'Tunu'}; % choose from 'Tunu', 'B40', 'NGRIP'
eruptions = {'UE 1453'};% list eruption(s) to consider
niters = 10000;
thresholdMIF = 0.1; % This is the number above which a sample is considered to have stratospheric sulfate
d34trange = 15; % range of d34S of tropospheric volcanic sulfate
d34tmin = -5; % minimum value of the range of d34S of volcanic sulfate
stratmin = 0; %starting minimum value of the d34S stratospheric
stratmax = 30; %starting maximum value of the d34S stratospheric
d34stratrestricted = true; % true if the strictest restrictions on d34S stratospheric are kept such that there is a monotonic decrease for d34Sstrat. This should be set to false for Huaynaputina

%% Run analyses core by core, and first determine what the background and uncertainty for the core should be

for kk = 1:length(cores)
    % find all data from a given core, and make sure it is sorted from
    % deepest to shallowest
    core_ind = find(strcmp(imported_data.Core(:,1),cores(kk)));
    core_data = sortrows(imported_data(core_ind,:), {'TopDepth'}, 'descend');

    %Determine the how variable the background is for the core by
    % taking standard deviation of concentration and d34S of all background
    % samples in the core
    if not(strcmp(cores(kk),'NGRIP'))
        all_bkgd_ind = find(strcmp(core_data.Type(:,1), 'bkgd')); %indices of all background samples
        bkgd_err = std(core_data.Conc(all_bkgd_ind)); % standard deviation of concentration of all background samples in core
        d34bkgd_err = std(core_data.d34S(all_bkgd_ind)); %standard deviation of d34S of all background samples in core
    else
        bkgd_err = 17; %standard deviation of NGRIP background samples run in STAiG lab
        d34bkgd_err = 1.3; %standard deviation of NGRIP background samples run in STAiG lab
    end

    for jj = 1:length(eruptions) % treat each eruption event separately
        % find all samples from a given event and save them in working matrix D
        eruption_ind = find(strcmp(core_data.Eruption(:,1),eruptions(jj)));
        D = core_data(eruption_ind, :);


        %find background sample(s)
        indbkgd = find(strcmpi(D.Type(:), 'bkgd'));

        if indbkgd % if there are background sample(s), calculate the mean of concentration and d34S
            bkgd =  mean(D.Conc(indbkgd)); %average background concentration for this eruption
            d34bkgd = mean(D.d34S(indbkgd)); % average d34S of background for this eruption
        elseif strcmp(cores(kk),'Tunu') && strcmp(eruptions(jj), 'UE 540') % for the case of 540 in Tunu there is no appropriate background so use the background from 536

            bkgd = 20.2; %concentration
            d34bkgd = 14.65; %isotope value
            disp('using Tunu 536 background for Tunu 540 event')

        else

            % if there is not a background
            disp('No background for:')
            disp(cores(kk))
            disp(eruptions(jj))
        end

        indevent = find(strcmpi(D.Type(:), 'event')); % Find all of the samples classified as having volcanic sulfate (i.e.'event')

        j = 1; %start at first event sample

        while j <= length(indevent)

            if D.D33S(indevent(j)) < thresholdMIF % if there is no measurable MIF signal

                %Just calculate the fvolc and its isotopic composition
                % load in the isotope data for this sample, and create Monte
                % Carlo matrices
                d34M = randn(niters,1)*D.d34S_err(indevent(j))/2+D.d34S(indevent(j));
                d33M = randn(niters,1)*D.d33S_err(indevent(j))/2+D.d33S(indevent(j));

                % Calculate the Monte Carlo matrices for the isotopic
                % composition of background sulfate
                d34b = randn(niters,1)*d34bkgd_err+d34bkgd;
                d33b = 1000 * ((1 + d34b./1000).^0.515 - 1);

                % calculate the fraction background and fraction volcanic and create Monte Carlo
                % matrices



                bkgdMC = randn(niters,1)*bkgd_err+bkgd;


                fb = bkgdMC./D.Conc(indevent(j)); %MC matrix for background fraction

                indexneg = find(fb(:,1)<0); %any negative fb limited to 0
                fb(indexneg) = 0;
                index = find(fb(:,1)>1); % find any impossible background fractions (greater than 1)

                fb(index) = 1; % must be 1



                fvolc = 1-fb; %calculate fraction volcanic

                % calculate isotopic composition of volcanic sulfate

                for i = 1:niters
                    if fvolc(i)~=0
                        d34volc(i) = (d34M(i) - fb(i).*d34b(i))./fvolc(i);
                        d33volc(i) = (d33M(i) - fb(i).*d33b(i))./fvolc(i);
                        D33volc(i) = (d33volc(i)/1000 - ((d34volc(i)./1000+1).^0.515 -1))*1000;
                    else

                        d34volc(i) = NaN;
                        d33volc(i) = NaN;
                        D33volc(i) = NaN;
                    end
                end

                % store median and standard deviation
                D.fvolc{indevent(j)} = median(fvolc,'omitnan');
                D.fvolc_err{indevent(j)} = std(fvolc,'omitnan');
                D.d34volc{indevent(j)} = median(d34volc,'omitnan');
                D.d34volc_err{indevent(j)} = std(d34volc,'omitnan');
                D.d33volc{indevent(j)} = median(d33volc,'omitnan');
                D.d33volc_err{indevent(j)} = std(d33volc,'omitnan');
                D.D33volc{indevent(j)} = median(D33volc,'omitnan');
                D.D33volc_err{indevent(j)} = std(D33volc,'omitnan');


                j = j+1; %go to next sample
            else %if there is measureable MIF signal
                firststrat_ind = indevent(j); %store this index as the index for the first sample with stratospheric sulfur

                % load in the isotope data for this sample, and create Monte
                % Carlo matrices
                d34M = randn(niters,1)*D.d34S_err(indevent(j))/2+D.d34S(indevent(j));
                d33M = randn(niters,1)*D.d33S_err(indevent(j))/2+D.d33S(indevent(j));




                % Calculate the Monte Carlo matrices for the isotopic
                % composition of background sulfate
                d34b = randn(niters,1)*d34bkgd_err+d34bkgd;
                d33b = 1000 * ((1 + d34b./1000).^0.515 - 1);


                %Calculate the Monte Carlo matrices for the isotopic
                %composition of tropospheric volcanic sulfate

                d34t = rand(niters,1)*d34trange+d34tmin;
                d33t = 1000 * ((1 + d34t./1000).^0.515 - 1);


                % calculate the fraction background and fraction volcanic and create Monte Carlo
                % matrices



                bkgdMC = randn(niters,1)*bkgd_err+bkgd;


                fb = bkgdMC./D.Conc(indevent(j)); %MC matrix for background fraction

                indexneg = find(fb(:,1)<0); %any negative fb limited to 0
                fb(indexneg) = 0;
                index = find(fb(:,1)>1); % find any impossible background fractions (greater than 1)

                fb(index) = 1; % must be 1



                fvolc = 1-fb; %calculate fraction volcanic

                % calculate isotopic composition of volcanic sulfate

                for i = 1:niters
                    if fvolc(i)~=0
                        d34volc(i) = (d34M(i) - fb(i).*d34b(i))./fvolc(i);
                        d33volc(i) = (d33M(i) - fb(i).*d33b(i))./fvolc(i);
                        D33volc(i) = (d33volc(i)/1000 - ((d34volc(i)./1000+1).^0.515 -1))*1000;
                    else

                        d34volc(i) = NaN;
                        d33volc(i) = NaN;
                        D33volc(i) = NaN;
                    end
                end

                % store median and standard deviation
                D.fvolc{indevent(j)} = median(fvolc,'omitnan');
                D.fvolc_err{indevent(j)} = std(fvolc,'omitnan');
                D.d34volc{indevent(j)} = median(d34volc,'omitnan');
                D.d34volc_err{indevent(j)} = std(d34volc,'omitnan');
                D.d33volc{indevent(j)} = median(d33volc,'omitnan');
                D.d33volc_err{indevent(j)} = std(d33volc,'omitnan');
                D.D33volc{indevent(j)} = median(D33volc,'omitnan');
                D.D33volc_err{indevent(j)} = std(D33volc,'omitnan');


                l = randn(niters,1)*0.006+0.608; %allow there to be uncertainty on lambda


                [solutions] = fstrat_MC(d34M, d33M,fb, d34b,d33b,  d34t, d33t,  l,stratmin, stratmax);
                D.solutions{indevent(j)} = solutions;
                D.fstrat{indevent(j)} = median(solutions(:,1));
                D.fstrat_err{indevent(j)} = std(solutions(:,1));
                D.d34Sstrat{indevent(j)} = median(solutions(:,2));
                D.d34Sstrat_err{indevent(j)} = std(solutions(:,2));
                j = length(indevent) +1; % get out of the while loop
            end
        end
        %Now go through the rest of the event samples
        for ii = firststrat_ind+1:indevent(end)
            MCsolutions = datasample(D.solutions{ii-1},niters);

            % load in the isotope data for this sample, and create Monte
            % Carlo matrices
            d34M = randn(niters,1)*D.d34S_err(ii)/2+D.d34S(ii);
            d33M = randn(niters,1)*D.d33S_err(ii)/2+D.d33S(ii);
            D33M = (d33M./1000 - ((d34M./1000+1).^0.515-1))*1000; %calculate the D33 of each MC iteration



            % Calculate the Monte Carlo matrices for the isotopic
            % composition of background sulfate
            d34b = randn(niters,1)*d34bkgd_err+d34bkgd;
            d33b = 1000 * ((1 + d34b./1000).^0.515 - 1);

            %Calculate the Monte Carlo matrices for the isotopic
            %composition of tropospheric volcanic sulfate

            %d34t = MCsolutions(:,3); %pass through d34S tropospheric options consistent with previous sample solutions
            d34t = rand(niters,1)*d34trange+d34tmin;
            d33t = 1000 * ((1 + d34t./1000).^0.515 - 1);

            % calculate the fraction background and fraction volcanic and create Monte Carlo
            % matrices



            bkgdMC = randn(niters,1)*bkgd_err+bkgd;


            fb = bkgdMC./D.Conc(ii); %MC matrix for background fraction

            indexneg = find(fb(:,1)<0); %any negative fb limited to 0
            fb(indexneg) = 0;
            index = find(fb(:,1)>1); % find any impossible background fractions (greater than 1)

            fb(index) = 1; % must be 1



            fvolc = 1-fb; %calculate fraction volcanic

            % calculate isotopic composition of volcanic sulfate

            for i = 1:niters
                if fvolc(i)~=0
                    d34volc(i) = (d34M(i) - fb(i).*d34b(i))./fvolc(i);
                    d33volc(i) = (d33M(i) - fb(i).*d33b(i))./fvolc(i);
                    D33volc(i) = (d33volc(i)/1000 - ((d34volc(i)./1000+1).^0.515 -1))*1000;
                else

                    d34volc(i) = NaN;
                    d33volc(i) = NaN;
                    D33volc(i) = NaN;
                end
            end

            % store median and standard deviation
            D.fvolc{ii} = median(fvolc,'omitnan');
            D.fvolc_err{ii} = std(fvolc,'omitnan');
            D.d34volc{ii} = median(d34volc,'omitnan');
            D.d34volc_err{ii} = std(d34volc,'omitnan');
            D.d33volc{ii} = median(d33volc,'omitnan');
            D.d33volc_err{ii} = std(d33volc,'omitnan');
            D.D33volc{ii} = median(D33volc,'omitnan');
            D.D33volc_err{ii} = std(D33volc,'omitnan');




            l = MCsolutions(:,4); %pass through lamda options consistent with previous sample solutions

            if ~strcmp(eruptions(jj), 'Huaynaputina') %for non-Huaynaputina samples
                stratmax_input = MCsolutions(:,2);%each subsequent sample must have a d34Sstrat lower than the previous sample
                stratmin_input = -60*ones(niters,1); % the minimum basic constraint on d34S strat is -60 per mil
                if d34stratrestricted %if want to include strictest restrictions on d34S strat
                    negindices = find(D33M <= 0); % indices if the measured D33S of the Monte Carlo simulation is negative
                    posindices = find(D33M >0); % % indices if the measured D33S of the Monte Carlo simulation is positive
                    stratmin_input(posindices,1) = MCsolutions(posindices,3); %if D33S is positive the minimum d34S can be is the starting value of d34S tropospheric
                    stratmax_input(negindices,1) = min(MCsolutions(negindices,2),MCsolutions(negindices,3)); %if D33S is negative the maximum d34S must be lower than both the starting tropospheric sulfate and the previous d34S stratospheric sample
                end


            else
                stratmax_input = stratmax; %for Huaynaputina, because two eruptions involved, we relax this constraint and keep the max to 30 per mil
                stratmin_input = -60;
            end




            [solutions] = fstrat_MC(d34M, d33M,fb, d34b,d33b,  d34t, d33t,  l,stratmin_input, stratmax_input);
            D.solutions{ii} = solutions;

            D.fstrat{ii} = median(solutions(:,1));
            D.fstrat_err{ii} = std(solutions(:,1));
            D.d34Sstrat{ii} = median(solutions(:,2));
            D.d34Sstrat_err{ii} = std(solutions(:,2));





        end
        eruptionsplit = split(eruptions{jj});
        if d34stratrestricted
            filename = strcat(cores(kk),'_',eruptionsplit(end),'_restricted');
        else
            filename = strcat(cores(kk),'_',eruptionsplit(end));
        end

        save(filename{1},'D')

    end

end

