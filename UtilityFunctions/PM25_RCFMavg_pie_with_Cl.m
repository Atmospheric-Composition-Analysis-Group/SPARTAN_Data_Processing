
function [Pie_made] = PM25_RCFMavg_pie_with_Cl(PM25_total, species_plot, Cl, K, Site_cities,loc)
% === comment out on Oct 30, 2024, Haihui [on-going] ===
% temporarily avoid reporting fraction of FITR OC in pie chars
% Instead, reporting only RM as the sum of current OC and RM

% % Don't want to include OC when only few filters has OC data
% ind = find(isnan(species_plot(:,end-1)) == 1 ); % OC is not available
% if length(ind) == size(species_plot,1)
%     species_plot(:,end-1) = 0; % no OC available, do not include in pie chart
% elseif length(ind) > 0.9*length(species_plot) && length(species_plot) < 100 % There are OC but only a few OC 
%     species_plot(:,end-1) = 0; % Don't include for pie chart
% else
%     species_plot(ind,end) = nan; % there are some OC data, mask out RM when OC not available to avoid non-representive OC-RM ratio
% end
% =============================================

% calculate
species = mean(species_plot ,1,'omitnan');
if isnan(species(end))
    species(end) = 0;
end

if sum(isnan(species)) == 0

    % SpecName = {'Water', 'SO4', 'NH4', 'NO3', 'Na', 'Dust', 'TEO', 'BC', 'OC', 'OM', 'Cl', 'K'}; % 11 x 1
    SpecName = {'Water', 'SO4', 'NH4', 'NO3', 'Na', 'Dust', 'TEO', 'BC',  'OM', 'Cl', 'K'}; % 10 x 1, no OC
    % --------------------
    Colors = [ 109,207,246; % water
        237,48,41;  % red sulfate
        240,103,166; % pink ammonium
        245,126,32; % orange nitrate
        57,84,165; % blue sodium
        252,238,30; % yellow dust
        128,130,133; % grey TEO
        35,31,32; % black carbon
        80, 184, 72; % green residual        55, 98, 60; % dark green OC
        56, 170, 165 ; % Cl
        180, 132, 212; % purple K
        ]./255;

    % change the order:
    species = [species, Cl, K];
    species = species(:,[8 7 6 5 10 11 4 3 2 1 9 ]);
    SpecName = SpecName([8 7 6 5 10 11 4 3 2 1 9 ]);
    Colors = Colors([8 7 6 5 10 11 4 3 2 1 9 ],:); % change the order

    % 2022-05-15 Haihui: replace negative species (usually OM) with 0
    notes ='';
    Ind = find(species<0);
    for i = 1:length(Ind) % prepare note in the chart
        notes = sprintf('%s\n%s = %.2f and is set to 0 in the chart.',notes,SpecName{Ind(i)},species(Ind(i)));
    end
    species(Ind)=0; % set negative to 0;
    % 20224-09-19 Haihui: replace nan with 0
    Ind = find(isnan(species)==1);
    for i = 1:length(Ind) % prepare note in the chart
        notes = sprintf('%s\n%s is nan and is set to 0 in the chart.', notes, SpecName{Ind(i)});
    end
    species(Ind) = 0; % set negative to 0;

    % Colors =
    % {'#A2142F','#0072BD','#D95319','k','#77AC30','#EDB120','#4DBEEE',''}; % AMS colors
    % SO4         NO3     NH4      BC    OM       Dust        SS
    figure(1)
    h = pie(species); % this prints fractions
    % h = pie(species,SpecName); % this prints species names
    htext = findobj(h,'Type','text');
    hPatch = findobj(h,'Type','patch');
    for sp = 1:length(SpecName)
        hPatch(sp).FaceColor = Colors(sp,:);
        htext(sp).String='';
        %     htext(sp).FontSize = 20;
        hold on
    end

    % centering the PM value in middle
    PM25str = sprintf('%1.0f',PM25_total);
    if PM25_total < 100
        fz = 44;
        sc = 0.11;
    else
        fz = 40;
        sc = 0.12;
    end
    plot(0, 0,'o', 'MarkerSize',90,'MarkerEdgeColor','k','MarkerFaceColor','w'); % PM2.5 white circle in centre
    text(-sc*length(PM25str),0,sprintf('%1.0f',PM25_total),'fontsize',fz,'fontweight','bold')
    hold on

    % 2022-05-15 Haihui: adding notes about negative number
    if ~isempty(notes)
        text(-0.1, 0, notes,'units','normalized')
    end
    % --------------------

    hold off
    T = title(sprintf('%s',Site_cities{loc}),'FontSize',30,'FontWeight','bold');
    set(T,'position',[-0.0039 1.05 1.001]);

    legend(SpecName, 'location','eastoutside')

    Pie_made = 1;
else
    Pie_made = 0;
end
