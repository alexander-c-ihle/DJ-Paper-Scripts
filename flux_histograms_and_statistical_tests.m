%% Create histogram figures of fluxes broken down by predicitive feature. 
%% Also does statistical tests to test relationship between fluxes and

clear
close all

% load data
data = readtable('all_instrument_compiled_fluxes_v2.xlsx');
BearCrk1 = readtable('BearCrk1_fluxes.xlsx');
CentrlPk = readtable('CentrlPk_fluxes.xlsx');
BearCrk2 = readtable('BearCrk2_fluxes.xlsx');
Chataqu2 = readtable('Chataqu2_fluxes.xlsx');
Chataq3 = readtable('Chataq3_fluxes.xlsx');
ClearSpr = readtable('ClearSpr_fluxes.xlsx');
Fthills1 = readtable('Fthills1_fluxes.xlsx');
Greeley1 = readtable('Greeley1_fluxes.xlsx');
Greeley2 = readtable('Greeley2_fluxes.xlsx');
Greeley3 = readtable('Greeley3_fluxes.xlsx');
MtHrmRd2 = readtable('MtHrmRd2_fluxes.xlsx');
PhilMill = readtable('PhilMill_fluxes.xlsx');
RedRckCan = readtable('RedRckCan_fluxes.xlsx');
RedRocks = readtable('RedRocks_fluxes.xlsx');
SKTrnst3 = readtable('SKTRnst3_fluxes.xlsx');
SKTrnst5 = readtable('SKTrnst5_fluxes.xlsx');
SKTrnst6 = readtable('SKTrnst6_fluxes.xlsx');
SodaLak3 = readtable('SodaLak3_fluxes.xlsx');
TomWatsn = readtable('TomWatsn_fluxes.xlsx');
TantraPk = readtable('TantraPk_fluxes.xlsx');

%% Seasonal comparason
x_winter = [];
x_summer = [];

% Seperate all winter vs. summer measurements
for i = 1: height(data)
    if data.winter(i) == 1;
    x_winter(i) = asinh(data.CH4_Flux__mg_CH4_m_2_day_1_(i));
    else x_winter(i) = "";
    end
   if data.winter(i) == 0;
      x_summer(i) = asinh(data.CH4_Flux__mg_CH4_m_2_day_1_(i));
   else x_summer(i) = "";
    end
end
x_summer = rmmissing(x_summer);
x_winter = rmmissing(x_winter);
x_summer_untransformed = sinh(x_summer);
x_winter_untransofmred = sinh(x_winter);

% Two-sided Wilcoxon rank test for statiscal significance of entire
% seasonal dataset
[wilcoxon_full_seasonal_pvalue, wilcoxon_full_seasonal_h, wilcoxon_full_seasonal_stats] = ranksum(x_summer, x_winter);

% Two-sided Wilcoxon rank test for only sites sampled in both summer and winter 
repeat_summer = [BearCrk1.summer_flux; CentrlPk.summer_flux; BearCrk2.summer_flux; Chataqu2.summer_flux; Chataq3.summer_flux; ...
    ClearSpr.summer_flux; Fthills1.summer_flux; Greeley1.summer_flux; Greeley2.summer_flux; Greeley3.summer_flux; MtHrmRd2.summer_flux; ...
    PhilMill.summer_flux; RedRckCan.summer_flux; RedRocks.summer_flux; SKTrnst3.summer_flux; SKTrnst5.summer_flux; SKTrnst6.summer_flux; SodaLak3.summer_flux; ...
    TomWatsn.summer_flux; TantraPk.summer_flux];
repeat_summer = rmmissing(repeat_summer);

repeat_winter = [BearCrk1.winter_flux; CentrlPk.winter_flux; BearCrk2.winter_flux; Chataqu2.winter_flux; Chataq3.winter_flux; ...
    ClearSpr.winter_flux; Fthills1.winter_flux; Greeley1.winter_flux; Greeley2.winter_flux; Greeley3.winter_flux; MtHrmRd2.winter_flux; ...
    PhilMill.winter_flux; RedRckCan.winter_flux; RedRocks.winter_flux; SKTrnst3.winter_flux; SKTrnst5.winter_flux; SKTrnst6.winter_flux; SodaLak3.winter_flux; ...
    TomWatsn.winter_flux; TantraPk.winter_flux];
repeat_winter = rmmissing(repeat_winter);

[wilcoxon_repeat_pvalue, wilcoxon_repeat_h, wilcoxon_repeat_stats] = ranksum(repeat_summer, repeat_winter);

% Two-sided Wilcoxon rank tests comparing individual site's summer and
% winter measurements
[wilcoxon_BearCrk1_pvalue, wilcoxon_BearCrk1_h, wilcoxon_BearCrk1_stats] = ranksum(BearCrk1.summer_flux, BearCrk1.winter_flux);
[wilcoxon_CentrlPk_pvalue, wilcoxon_CentrlPk_h, wilcoxon_CentrlPk_stats] = ranksum(CentrlPk.summer_flux, CentrlPk.winter_flux);
[wilcoxon_BearCrk2_pvalue, wilcoxon_BearCrk2_h, wilcoxon_BearCrk2_stats] = ranksum(BearCrk2.summer_flux, BearCrk2.winter_flux);
[wilcoxon_Chataqu2_pvalue, wilcoxon_Chataqu2_h, wilcoxon_Chataqu2_stats] = ranksum(Chataqu2.summer_flux, Chataqu2.winter_flux);
[wilcoxon_Chataq3_pvalue, wilcoxon_Chataq3_h, wilcoxon_Chataq3_stats] = ranksum(Chataq3.summer_flux, Chataq3.winter_flux);
[wilcoxon_ClearSpr_pvalue, wilcoxon_ClearSpr_h, wilcoxon_ClearSpr_stats] = ranksum(ClearSpr.summer_flux, ClearSpr.winter_flux);
[wilcoxon_Fthills1_pvalue, wilcoxon_Fthills1_h, wilcoxon_Fthills1_stats] = ranksum(Fthills1.summer_flux, Fthills1.winter_flux);
[wilcoxon_Greeley1_pvalue, wilcoxon_Greeley1_h, wilcoxon_Greeley1_stats] = ranksum(Greeley1.summer_flux, Greeley1.winter_flux);
[wilcoxon_Greeley2_pvalue, wilcoxon_Greeley2_h, wilcoxon_Greeley2_stats] = ranksum(Greeley2.summer_flux, Greeley2.winter_flux);
[wilcoxon_Greeley3_pvalue, wilcoxon_Greeley3_h, wilcoxon_Greeley3_stats] = ranksum(Greeley3.summer_flux, Greeley3.winter_flux);
[wilcoxon_MtHrmRd2_pvalue, wilcoxon_MtHrmRd2_h, wilcoxon_MtHrmRd2_stats] = ranksum(MtHrmRd2.summer_flux, MtHrmRd2.winter_flux);
[wilcoxon_PhilMill_pvalue, wilcoxon_PhilMill_h, wilcoxon_PhilMill_stats] = ranksum(PhilMill.summer_flux, PhilMill.winter_flux);
[wilcoxon_RedRckCan_pvalue, wilcoxon_RedRckCan_h, wilcoxon_RedRckCan_stats] = ranksum(RedRckCan.summer_flux, RedRckCan.winter_flux);
[wilcoxon_RedRocks_pvalue, wilcoxon_RedRocks_h, wilcoxon_RedRocks_stats] = ranksum(RedRocks.summer_flux, RedRocks.winter_flux);
[wilcoxon_SKTrnst3_pvalue, wilcoxon_SKTrnst3_h, wilcoxon_SKTrnst3_stats] = ranksum(SKTrnst3.summer_flux, SKTrnst3.winter_flux);
[wilcoxon_SKTrnst5_pvalue, wilcoxon_SKTrnst5_h, wilcoxon_SKTrnst5_stats] = ranksum(SKTrnst5.summer_flux, SKTrnst5.winter_flux);
[wilcoxon_SKTrnst6_pvalue, wilcoxon_SKTrnst6_h, wilcoxon_SKTrnst6_stats] = ranksum(SKTrnst6.summer_flux, SKTrnst6.winter_flux);
[wilcoxon_SodaLak3_pvalue, wilcoxon_SodaLak3_h, wilcoxon_SodaLak3_stats] = ranksum(SodaLak3.summer_flux, SodaLak3.winter_flux);
[wilcoxon_TomWatsn_pvalue, wilcoxon_TomWatsn_h, wilcoxon_TomWatsn_stats] = ranksum(TomWatsn.summer_flux, TomWatsn.winter_flux);
[wilcoxon_TantraPk_pvalue, wilcoxon_TantraPk_h, wilcoxon_TantraPk_stats] = ranksum(TantraPk.summer_flux, TantraPk.winter_flux);

% Construct table compiling stats for every site sampled in both summer and
% winter
site_name = ["all repeat summer/winter sites";"BearCrk1";"CentrlPk"; "BearCrk2";"Chataqu2"; "Chataq3"; "ClearSpr";"Fthills1";...
    "Greeley1";"Greeley2";"Greeley3";"MtHrmRd2";"PhilMill"; "RedRckCan";"RedRocks";"SKTrnst3";"SKTrnst5";"SKTrnst6";"SodaLak3";"TomWatsn";"TantraPk"];
pvalue=  [wilcoxon_repeat_pvalue; wilcoxon_BearCrk1_pvalue; wilcoxon_CentrlPk_pvalue;wilcoxon_BearCrk2_pvalue;wilcoxon_Chataqu2_pvalue; wilcoxon_Chataq3_pvalue;wilcoxon_ClearSpr_pvalue; wilcoxon_Fthills1_pvalue; wilcoxon_Greeley1_pvalue;wilcoxon_Greeley2_pvalue;wilcoxon_Greeley3_pvalue; wilcoxon_MtHrmRd2_pvalue; wilcoxon_PhilMill_pvalue;wilcoxon_RedRckCan_pvalue;wilcoxon_RedRocks_pvalue;wilcoxon_SKTrnst3_pvalue;wilcoxon_SKTrnst5_pvalue;wilcoxon_SKTrnst6_pvalue;wilcoxon_SodaLak3_pvalue;wilcoxon_TomWatsn_pvalue; wilcoxon_TantraPk_pvalue];
h = [wilcoxon_repeat_h; wilcoxon_BearCrk1_h; wilcoxon_CentrlPk_h;wilcoxon_BearCrk2_h; wilcoxon_Chataqu2_h;wilcoxon_Chataq3_h;wilcoxon_ClearSpr_h; wilcoxon_Fthills1_h;wilcoxon_Greeley1_h;wilcoxon_Greeley2_h; wilcoxon_Greeley3_h; wilcoxon_MtHrmRd2_h;wilcoxon_PhilMill_h;wilcoxon_RedRckCan_h;wilcoxon_RedRocks_h;wilcoxon_SKTrnst3_h;wilcoxon_SKTrnst5_h; wilcoxon_SKTrnst6_h;wilcoxon_SodaLak3_h; wilcoxon_TomWatsn_h; wilcoxon_TantraPk_h];
seasonal_comparason = table(site_name, pvalue, h);

% Plot summer vs. winter fluxes
figure('color','white')
subplot(2,2,1);
set(subplot(2,2,1),'Color','White')
hold on;
histogram(x_summer,[-3:1/3:6],'FaceColor',"#5C7F1C",EdgeColor="#5e5e5e", FaceAlpha = 0.8);
yscale log
ylim ([0 400]);
yticklabels({'1', '10', '100'});
histogram(x_winter,[-3:1/3:6],'FaceColor',"#23484F",EdgeColor="#5e5e5e", FaceAlpha = 1);
ax1 = gca;
ax1.Color = 'white';
ax1.XColor = 'black';
ax1.YColor = 'black';
ax1.FontSize = 24;
ax1.FontName = 'calibri';
lgd1 = legend('Summer', 'Winter', 'Color', 'white', 'TextColor', 'black');
lgd1.EdgeColor = 'black';
xlabel('Flux (mg CH_4 m^{-2} day^{-1})');ylabel('Measurement count');
xticks([-2.99822295,  -0.881373587, 0, 0.881373587, 2.99822295,    5.2983 ]) % Asinh of x axis tick label values
xticklabels({'-10','-1','0','1','10','100'})
text(ax1, -0.2, 1.25, 'a)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'calibri', 'FontWeight', 'bold','Color', 'black', 'VerticalAlignment', 'top');
hold off;


%% E vs. W of fold axis comparason
w_of_axis = [];
e_of_axis = [];
dummy = []; % Dummy variable for plotting overlap of histograms on legend

% Seperate E and W flux measurements
for i = 1:height(data)
      if data.w_of_fold_axis(i) == 1
          w_of_axis(i) = asinh(data.CH4_Flux__mg_CH4_m_2_day_1_(i));
      else w_of_axis(i) = "";
      end
     if data.w_of_fold_axis(i) == 0
          e_of_axis(i) = asinh(data.CH4_Flux__mg_CH4_m_2_day_1_(i));
      else e_of_axis(i) = "";
     end 
end 


% Remove NaNs
w_of_axis = rmmissing(w_of_axis);
e_of_axis = rmmissing(e_of_axis);

% Make untransformed versions of e, w data
w_of_axis_untransformed = sinh(w_of_axis);
e_of_axis_untransformed = sinh(e_of_axis);

% Subset of E vs. W of fold axis data with only positive fluxes above
% detection limit
e_of_axis_allpos = [];
for i = 1:length(e_of_axis_untransformed)
      if e_of_axis_untransformed(i) >=0.05
          e_of_axis_allpos(i) = (e_of_axis_untransformed(i));
        else 
        e_of_axis_allpos(i) = "NaN";
      end     
end 

w_of_axis_allpos = [];
for i = 1:length(w_of_axis_untransformed)
      if w_of_axis_untransformed(i) >=0.05
          w_of_axis_allpos(i) = (w_of_axis_untransformed(i));
        else 
        w_of_axis_allpos(i) = "NaN";
      end     
end

% Remove NaNs, calculate stats
w_of_axis_allpos = rmmissing(w_of_axis_allpos);
e_of_axis_allpos = rmmissing(e_of_axis_allpos);
w_of_axis_allpos_mean = mean(w_of_axis_allpos);
w_of_axis_allpos_median = median(w_of_axis_allpos);
w_of_axis_allpos_standard_error = std(w_of_axis_allpos)/sqrt(width(w_of_axis_allpos));
e_of_axis_allpos_mean = mean(e_of_axis_allpos);
e_of_axis_allpos_median = median(e_of_axis_allpos);
e_of_axis_allpos_standard_error = std(e_of_axis_allpos)/sqrt(width(e_of_axis_allpos));

% Two-sided Wilcoxon rank test for statiscal significance
[wilcoxon_all_e_w_pvalue, wilcoxon_all_e_w_h, wilcoxon_all_e_w_stats] = ranksum(w_of_axis,e_of_axis);
[wilcoxon_allpos_pvalue, wilcoxon_allpos_h, wilcoxon_allpos_stats] = ranksum(e_of_axis_allpos,w_of_axis_allpos);

% Plot
subplot(2,2,2);
set(subplot(2,2,2),'Color','White')
hold on
histogram(e_of_axis,[-3:1/3:6], FaceColor = "#386E38",EdgeColor="#5e5e5e", FaceAlpha = 0.7);%AF7443
histogram(w_of_axis,[-3:1/3:6],FaceColor = "#CFEAEF", EdgeColor="#5e5e5e", FaceAlpha = 0.7);
histogram(dummy, FaceColor = "#B4D2CA", EdgeColor="#5e5e5e");
yscale log
ylim ([0 400]);
yticklabels({'1', '10', '100'});
ax2 = gca;
ax2.Color = 'white';
ax2.XColor = 'black';
ax2.YColor = 'black';
ax2.FontSize = 24;
ax2.FontName = 'calibri';
lgd2 = legend('E of axis', 'W of axis', 'E/W of axis overlap', 'Color', 'white', 'TextColor', 'black');
lgd2.EdgeColor = 'black';
xlabel('Flux (mg CH_4 m^{-2} day^{-1})');ylabel('Measurement count');%title('inverse hyperbolic sin by season');
legend();
xticks([-2.99822295,  -0.881373587, 0, 0.881373587, 2.99822295,    5.2983 ])
xticklabels({'-10','-1','0','1','10','100'})
text(ax2, -0.2, 1.25, 'b)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'calibri', 'FontWeight', 'bold','Color', 'black', 'VerticalAlignment', 'top');
hold off;


%% On vs. off geologic units of interest comparason
on_geo = [];
off_geo = [];

for i = 1:height(data)
      if data.on_geo_of_interest(i) == 1
          on_geo(i) = asinh(data.CH4_Flux__mg_CH4_m_2_day_1_(i));
        else 
        on_geo(i) = "";
      end
      if data.on_geo_of_interest(i) == 0
          off_geo(i) = asinh(data.CH4_Flux__mg_CH4_m_2_day_1_(i));
        else 
        off_geo(i) = "";
      end     
end 

% Remove NaNs
on_geo = rmmissing(on_geo);
off_geo = rmmissing(off_geo);

% Make untransformed versions of data
on_geo_untransformed = sinh(on_geo);
off_geo_untransformed = sinh(off_geo);

% Subset of data with only positive fluxes above detection limit
on_geo_allpos = [];
for i = 1:length(on_geo_untransformed)
      if on_geo_untransformed(i) >=0.05
          on_geo_allpos(i) = (on_geo_untransformed(i));
        else 
        on_geo_allpos(i) = "NaN";
      end     
end 

off_geo_allpos = [];
for i = 1:length(off_geo_untransformed)
      if off_geo_untransformed(i) >=0.05
          off_geo_allpos(i) = (off_geo_untransformed(i));
        else 
        off_geo_allpos(i) = "NaN";
      end     
end

% remove NaNs, calculate stats
on_geo_allpos = rmmissing(on_geo_allpos);
off_geo_allpos = rmmissing(off_geo_allpos);
on_geo_allpos_mean = mean(on_geo_allpos);
on_geo_allpos_median = median(on_geo_allpos);
on_geo_allpos_standard_error = std(on_geo_allpos)/sqrt(width(on_geo_allpos));
off_geo_allpos_mean = mean(off_geo_allpos);
off_geo_allpos_median = median(off_geo_allpos);
off_geo_allpos_standard_error = std(off_geo_allpos)/sqrt(width(off_geo_allpos));

% Do two-sided Wilcoxon rank test for statiscal significance
[wilcoxon_on_off_geo_pvalue, wilcoxon_on_off_geo_h, wilcoxon_on_off_geo_stats] = ranksum(off_geo_untransformed,on_geo_untransformed);
[wilcoxon_allpos__on_off_geo_pvalue, wilcoxon_allpos_on_off_geo_h, wilcoxon_allpos_on_off_geo_stats] = ranksum(off_geo_allpos,on_geo_allpos);

% Plot
subplot(2,2,3);
hold on
histogram(off_geo,[-3:1/3:6],FaceColor = "#C8A5F7",EdgeColor="#5e5e5e", FaceAlpha = 1);
histogram(on_geo,[-3:1/3:6],FaceColor = "#201F38", EdgeColor="#5e5e5e" , FaceAlpha = 1);
yscale log
ylim ([0 400]);
yticklabels({'1', '10', '100'});
ax3 = gca;
ax3.Color = 'white';
ax3.XColor = 'black';
ax3.YColor = 'black';
ax3.FontSize = 24;
ax3.FontName = 'calibri';
lgd3 = legend('Off hydrocarbon unit(s)', 'On hydrocarbon unit(s)', 'Color', 'white', 'TextColor', 'black');
lgd3.EdgeColor = 'black';
xlabel('Flux (mg CH_4 m^{-2} day^{-1})');ylabel('Measurement count');
xticks([-2.99822295,  -0.881373587, 0, 0.881373587, 2.99822295,    5.2983 ])
xticklabels({'-10','-1','0','1','10','100'})
text(ax3, -0.2, 1.25, 'c)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'calibri', 'FontWeight', 'bold','Color', 'black', 'VerticalAlignment', 'top');
hold off;

%% Flux magnitude vs distance to W basin boundary
% Do Spearman Rank test for statistical significance
[rho_all, pval_all] = corr(data.distance_to_w_boundary_m, data.CH4_Flux__mg_CH4_m_2_day_1_, 'Type', 'Spearman');

% do same test on just positive fluxes
both_fluxes_allpos = [];
w_boundary_distance_allpos = [];

for i = 1:length(data.CH4_Flux__mg_CH4_m_2_day_1_)
    if data.CH4_Flux__mg_CH4_m_2_day_1_(i) >= 0.05
        both_fluxes_allpos(end+1) = data.CH4_Flux__mg_CH4_m_2_day_1_(i); % Append positive flux
        w_boundary_distance_allpos(end+1) = data.distance_to_w_boundary_m(i); % Append corresponding distance
    end
end

% Convert to column vectors 
both_fluxes_allpos = both_fluxes_allpos(:);
w_boundary_distance_allpos = w_boundary_distance_allpos(:);

% Test positive fluxes
[rho_allpos, pval_allpos] = corr(w_boundary_distance_allpos, both_fluxes_allpos, 'Type', 'Spearman');
subplot(2,2,4);
scatter1 = scatter(data.distance_to_w_boundary_m/1000,asinh(data.CH4_Flux__mg_CH4_m_2_day_1_),'MarkerEdgeColor','#000000'); %#B9786C

% Plot
ax4 = gca;
ax4.Color = 'white';
ax4.XColor = 'black';
ax4.YColor = 'black';
ax4.FontSize = 24;
ax4.FontName = 'calibri';
yline(0, '--r', 'LineWidth', 2);
lgd4 = legend('', 'y = 0', 'Color', 'white', 'TextColor', 'black');
lgd4.EdgeColor = 'black';
xlabel('Distance to W basin boundary (km)');ylabel('Flux (mg CH_4 m^{-2} day^{-1})');%title('inverse hyperbolic sin by season');


yticks([-2.99822295,   0,  2.99822295,    5.2983 ])
yticklabels({'-10','0','10','100'})
text(ax4, -0.2, 1.25, 'd)', 'Units', 'normalized', 'FontSize', 24, 'FontName', 'Calibri', 'FontWeight', 'bold','Color', 'black', 'VerticalAlignment', 'top');
hold off;



