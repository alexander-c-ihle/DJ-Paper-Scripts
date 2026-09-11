%% sample synthetic datasets to determine flux limits
close all
clear all
addpath '.........\utils'

%% parameters

% grid limits and inncrements (degrees)
xlims = [-105.5 -102.5];
ylims = [38 41];

% spatial correlation/smoothing length scale
s_scale = 10000; % in m.

% Rejection tolerance. Fixes edge problem. Any of the simulations whose random 0-1 numbers 
% are >0.995 AND lie outside the basin area will be thrown out.
peak_tol = 0.995; 

% flux parameters
Freg_range = [0 10];

% looping
nv = [1:5]; % number of versions (use when the code takes a long time)
nsim = 10; % number of monte carlo interations per version

% || Basin Specific||
% Conversion factor for grid resolution - linearly scales with SC
% to maintain 3 grid cells per hotspot at 100 m SC and 10 grid cells per
% hotspot at 100 km SC. For SC lenghts between these two extremes, the number
% of grid cells per hotspot feature linearly scales with SC in the
% following four lines.
xinc1 = 3e-4; xinc2 = 0.09; % xinc, grid resolution for that SC. 
SC1 = 100; SC2 = 100000; % spatial correlation scale (SC)
CFm = ((xinc2-xinc1)/(SC2-SC1)); % conversion factor (SLOPE Y= MX + B) from SC to number of boxes
CFb = xinc1 - CFm*SC1;

%% load the flux data

% read the table fo flux data
T = readtable(['/all_instrument_compiled_fluxes_v2.csv']);

% unpack data, fluxes in mg ch4 m^-2 day^-1
samples.Fch4 = T.("CH4_Flux__mg_CH4_m_2_day_1_");samples.Fch4_err = T.("CH4_Flux_Uncertainty__mg_CH4_m_2_day_1_");
samples.bio_prob = T.("bio_prob");
samples.lon = T.("Latitude__DD_");samples.lat = T.("Longitude__DD_");
samples.date = T.("Date");

N = find(samples.Fch4<0); % Index negative fluxes for exclusion

% load in basin boundary
S = shaperead('DJ_clipped_shapefile_Projected.shp','UseGeoCoords',true); % Read in the shapefile.
lat = S.Lat; lon = S.Lon; % Decimal Degrees, Define Latitude and Longitude from the structure create from the shapefile's attributes.

%% set up grid

% set grid resolution
xinc = s_scale*CFm + CFb;
yinc = s_scale*CFm + CFb;

% gridcell boundaries
xbnd = xlims(1):xinc:xlims(2);
ybnd = ylims(1):yinc:ylims(2);

% gridcell centers
xcen = (xbnd(1:end-1)+xbnd(2:end))/2;
ycen = (ybnd(1:end-1)+ybnd(2:end))/2;
[X,Y] = meshgrid(xcen,ycen);

% gridcell area
[area,dx,dy] = grid_area(xcen,ycen,xinc,yinc);
[ny,nx] = size(area);

% index grid cells within basin
[in,out] = inpolygon(X,Y,lon,lat);
area(~in) = NaN;

% define smoothing windows. dx is distance between centers on longitude,
% which changes, dy is distance between centers on latitude, which doesn't
% change
xwin = round(s_scale/min(dx(:)));
ywin = round(s_scale/dy);

% Minimum flux fraction (given size of feature). The smallest size a
% hotspot can be divided by the total area. For low-flux scenarios you're
% testing, the code here scales down the intensity of hotspots so you can
% still test a reagoinal mean flux with at least one hotspot present.
fpos_min = pi*((.5*s_scale)^2)/sum(area(:),"omitmissing");

%% set up probability distributions

% Number of samples. May have to increase/decrease to populate entire study
% area
nsamp = 1e6;

% load Etiope positive distribution
load etiope_pos_pdf.mat pdfun F_etio
Fmax = max(F_etio);
pd_etio = pdfun;
fpool = sinh(random(pd_etio,nsamp,1));
fpool = fpool(fpool<Fmax);
pd_etio_mu = mean(fpool);

%% sample synthetic datasets

% Versions to save out (optional)
for v = min(nv):max(nv)
    for i = 1:nsim
        if mod(i,10)==0; disp(['iteration ' num2str(i)]); end

        % step 1: create a valid noise map
        valmap = 0;
        while valmap==0
            % generate padded white noise
            rn = rand(ny+2*ywin,nx+2*xwin);
            % smooth and scale the noise
            smn_tmp = smoothxy2(rn,xwin,ywin);
            smn_range = max(smn_tmp(:)) - min(smn_tmp(:));
            smn_tmp = (smn_tmp - min(smn_tmp(:)))/smn_range;
            % cut to size    
            smn = smn_tmp((ywin+1):(end-ywin),(xwin+1):(end-xwin));
            smn(~in) = 0; % cut to basin boundary
            % is map valid?
            if squeeze(max(max(smn,[],1),[],2))>peak_tol
                valmap=1;
            end
        end
        [jnk,Is] = sort(smn(:));
        Ival_jnk = find(~isnan(jnk));
        Ival = Is(Ival_jnk);
    
        % step 2: get target regional mean from uniform distribution
        Freg_target(i) = min(Freg_range) + rand(1,1)*(max(Freg_range)-min(Freg_range));
        fpos = max(fpos_min,Freg_target(i)/pd_etio_mu);
        npos = round(fpos*length(Ival));
        
        % step 3: fill the noise map with fluxes
        F_synth = smn*NaN;
        randraw = randsample(fpool,npos);
        randraw = [zeros(length(Ival)-npos,1);randraw];
        F_synth(Ival) = sort(randraw); % synthetic distribution
        
        % step 4: scale
        I = ~isnan(F_synth(:)); % Land values (takes out NaNs)
        F_synth_vert = F_synth(:); area_vert = area(:);
        F_synth = Freg_target(i)*F_synth/(F_synth(:)'*area(:)/sum(area(:)));
        area(~in) = 0;
        Freg(i) = (F_synth(:)'*area(:))/sum(area(:),"omitmissing"); % Accounts for differences
        
        % step 5: define sample set
        F_samp = samples.Fch4 + normrnd(0,1,length(samples.Fch4)).*samples.Fch4_err;
        Igeo = rand(length(samples.Fch4),1)>(samples.bio_prob/100);
        for j = 1:length(F_samp) % cut off uncertainty at twice the magnitude of the flux to remove bias from setting pdf below 0 to 0.
            if F_samp(j)>2*samples.Fch4(j);F_samp(j) = 2*samples.Fch4(j);else;end
        end
        F_samp(N) = 0; % Set negative fluxes to 0. (some samples, unc will make it positive)
        F_samp = max(0,F_samp(Igeo)); % Only positive, geologic fluxes are kept
    
        % step 6: calculate stats of sample set
        F_samp_mu(i) = mean(F_samp);
        F_samp_sd(i) = std(F_samp);
        F_se = std(F_samp)/sqrt(length(F_samp));
        F_upper(i) = F_samp_mu(i)+2*F_se;
        F_lower(i) = F_samp_mu(i)-2*F_se;
    
        % step 7: sample the synthetic distribution
        F_synth_samp = interp2(X,Y,F_synth,samples.lon(Igeo),samples.lat(Igeo));
        F_synth_mu(i) = mean(F_synth_samp,"omitnan");
        F_synth_sd(i) = std(F_synth_samp);
    end
  
    %% make output file    
    output.Freg = Freg;
    output.F_synth_mu = F_synth_mu;
    output.F_synth_sd = F_synth_sd;
    output.F_samp_mu = F_samp_mu;
    output.F_samp_sd = F_samp_sd;
    output.F_upper = F_upper;
    output.F_lower = F_lower;
   filename = sprintf("............/DJ_output"+"_SC"+string(s_scale)+"_v"+string(v)+".mat");
   save(filename, "output");
end

