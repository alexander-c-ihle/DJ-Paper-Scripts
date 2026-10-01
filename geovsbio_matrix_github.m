%% ------------ Geologic vs. Biologic Decision Tree -----------------

% This code is intended to interpret whether positive fluxes are
% thermogenic or biogenic in origin. The input dataset must be a finalized
% set of fluxes with full auxilary data. The output will contain the set of
% positive flux measurements with assigned probablity of being biogenic in
% origin. 

% It should be noted that definitive methods for determining source include
% geochemical (such as Carbon 13 and 14) measurements and analyses. These  
% are absent in these measurements and sites. This includes measurements of
% ethane along with methane from the Picarro analyzer. 
 

% Characteristics of the DJ:
%   - Deep gas production, so surface fluxes are not widespread and are generally small
%          
%   - Basin is dry, resulting in low soil moistures. Likely to have
%       few biologic fluxes 

% Housekeeping
clear
close all

%% Designate method matrix will assign % biologic likelihood
% option 1 with a linear interpolation
option = "1: Warm and Wet"; % "1: Warm and Wet", "2: Warm, Wet, and Glacial Till"
interp_meth = "linear"; % "linear" or "cubic"

%% Load Flux Data
T = readtable("all_instrument_compiled_fluxes_v2.csv");%,"Sheet","Combined_RAW");
f = T.CH4_Flux__mg_CH4_m_2_day_1_;

I1 = find(f>0.05); % Positive fluxes above our detection limit
I = sort(unique(I1));

data = T(I,:);
% Create seperate variables of all to columns of T that will be used in the
% decision tree/matrix 
smois = data.Soil_Moisture__VWC_;
% Change missing soil moisture values to 10%. Will effectivley designate
% positive fluxes without recorded soil moistures as geologic.
smois(isnan(smois)) = 10;
stemp = data.Soil_Temperature__C_;
bigwinter = data.bigger_winter_flux;
smallwinter = data.smaller_winter_flux;


geobio_categ = repmat("unknown", 1, length(data.Unique_Sample_ID));
bio_prob = (zeros(1,length(data.Unique_Sample_ID))); 

%% Create Probability Matrix
if option == "1: Warm and Wet"
    if interp_meth == "linear"
        [X,Y] = meshgrid(0:1);
        V = [10, 30; 50, 90];
        
        % Create query grid
        [Xq,Yq] = meshgrid(0:0.01:1);
        Vq = interp2(X,Y,V,Xq,Yq,"linear");
        
        figure(1)
        set(gcf,'Color','w');
        surf(Xq,Yq,Vq);
        title('Option 1: Linear, Warm & Wet');

    elseif interp_meth == "cubic"
        % Data/Original
        [X,Y] = meshgrid(0:2);
        V = [10, 20, 30; 50, 60, 75; 70, 85, 90];
        
        % Create query grid
        [Xq,Yq] = meshgrid(0:0.02:2);
        Vq = interp2(X,Y,V,Xq,Yq,"cubic");
        I = find(Vq>90);
        Vq(I) = 90; % Any value over 90 is set to 90

        figure(1)
        set(gcf,'Color','w');
        surf(Xq,Yq,Vq);
        title('Option 1: Cubic, Warm & Wet');
    end
    zlabel("Probability the Flux is Biological (%)")
    xlabel("Temperature")
    ylabel("Moisture")
end

%% Decision Tree
% A logical descision tree to make binary decisions on geologic/biologic
% (either/or). Remaining samples will be fit onto a probability matrix of
% being modern, biogenic fluxes.

x_cord_temp = (zeros(1,length(data.Unique_Sample_ID)))*NaN; 
y_cord_mois = (zeros(1,length(data.Unique_Sample_ID)))*NaN;
z_cord_geo = (zeros(1,length(data.Unique_Sample_ID)))*NaN;
var3 = (zeros(1,length(data.Unique_Sample_ID)))*NaN; % for plotting flux values on final matrix plot

for i = 1:length(data.Unique_Sample_ID)

    % Assign geologic flux to sites where fluxes were larger in winter
    % than in summer.
        if geobio_categ(i) == "unknown" && bigwinter(i) == 1
           geobio_categ(i) = "Geologic";
        end

    % Assign geologic fluxes to remaining uncategorized fluxes with soil
    % moistures less than 20%. 
        if geobio_categ(i) == "unknown" && smois(i) <= 20
           geobio_categ(i) = "Geologic";
        end
    % Assign geologic fluxes to remaining uncategorized fluxes where fluxes
    % in winter were smaller than in summer. The remaining uncategorized
    % fluxes will be sent to the matrix. 
        if geobio_categ(i) == "unknown" && smallwinter(i) == 1
           geobio_categ(i) = "Biologic";
           bio_prob(i) = 100;
        end
end

%% Apply matrix to any positive fluxes that made it through the decision
%% tree without being categorized as geologic
for i = 1:length(data.Unique_Sample_ID)
    if geobio_categ(i)=="unknown"
        if stemp(i) > 20
           stemp(i) = 20;
        end
        x_cord_temp(i) = round((stemp(i)/20)*length(Vq));
        y_cord_mois(i) = round((smois(i)/50)*length(Vq));
        var3(i) = data.CH4_Flux__mg_CH4_m_2_day_1_(i);
        % Coordinates can't be 0 in MATLAB
        if x_cord_temp(i) <= 0; x_cord_temp(i) = 1; end
        if y_cord_mois(i) == 0; y_cord_mois(i) = 1; end
        if option == "1: Warm and Wet"
            disp([i, x_cord_temp(i), y_cord_mois(i)])
            bio_prob(i) = Vq(x_cord_temp(i),y_cord_mois(i));
            bio_prob(i) = Vq(x_cord_temp(i),y_cord_mois(i));
        end
    end
end

point_color = (zeros(1,length(var3)))*NaN;
for i = 1:length(var3)
    if var3(i) > 1
        point_color(i) = 1; 
    else point_color(i) = 0;
    end 
end

% Plot scatter points with colors based on var3

% Define the two hexadecimal colors
color1 = '#FD8D3C';
color2 = '#FBCD5C';


% Convert hexadecimal to RGB
rgbColor1 = sscanf(color1(2:end), '%2x%2x%2x', [1 3]) / 255;
rgbColor2 = sscanf(color2(2:end), '%2x%2x%2x', [1 3]) / 255;

% Initialize color matrix
colors = zeros(length(point_color), 3);

% Assign colors based on the value of point_color
colors(point_color == 1, :) = repmat(rgbColor1, sum(point_color == 1), 1);
colors(point_color == 0, :) = repmat(rgbColor2, sum(point_color == 0), 1);


if option == "1: Warm and Wet"
    figure (2)
    set(gcf, 'Color', 'w');
    hold on
    % plot background color ramp
    h = pcolor(Vq);
    colormap('bone')
    a = colorbar;
    ylabel(a, '% Likelihood biologic flux');
    a.Color = 'black';          % Sets tick marks, tick labels, and outline to black
    a.Label.Color = 'black';    % Sets the ylabel text to black
    set(h, 'EdgeColor', 'none');% Remove grid lines 
    a.FontSize = 18;
    a.FontName = 'airial';
   
    % Exclude the pcolor plot from the legend
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
   
    % Plot scatter points
    s= scatter(x_cord_temp, y_cord_mois, 250, colors, 'square', 'filled');
    set(get(get(s, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    hold on
scatter(NaN, NaN, 150, rgbColor1, 'square', 'filled', ...
    'DisplayName', 'Flux ≥ 1 & < 10 mg CH_4 m^{-2} day^{-1}');
     % Add the legend, change marker size
  lgnd = legend('show','Location','best');
  lgnd.Color = 'white';        % Legend background color
  lgnd.TextColor = 'black';    % Legend text color
  lgnd.FontSize = 12;          % Font size
  lgnd.FontName = 'arial';     % Font name
     xlim([-2,102]);
    ylim([0,length(Vq)]);
    yticks(0:20:100);
    yticklabels(["","10","20","30","40","50"]);
    xticks(0:25:100);
    xticklabels(["0", "5","10","15","20"]);
    xlabel("Soil temperature (°C)");
    ylabel("Soil moisture (% VWC)");
    %colormap turbo
    ax = gca;
    ax.Color = 'white';
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.ZColor = 'black';
    ax.FontSize = 18;
    ax.FontName = 'arial';
    ax.TickDir = 'out';
    ax.XAxis.Color = 'black';
    ax.YAxis.Color = 'black';

else
end

