% This file is named: qdot_noise_display.m
% Displays maps for calculations in qdot_noise.m
% Saves the input qdot data and the resulted figure.
% Alternately, load data from qdot_noise.mat and run this.

% This script uses the bluewhitered color map from: www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
fighandle(1) = figure('rend','painters','pos',[900 500 1000 450]);

spdIdV = subplot(1,5,1);

dIdV = zeros(size(current));
dV = biasVv(2)-biasVv(1);

for Ei = 1:length(gateVv)
    for Vi = 2:length(biasVv)
        dIdV(Ei,Vi)= (current(Ei,Vi)-current(Ei,Vi-1))./dV; 
    end
end
ih=imagesc(gateVv,biasVv,real(dIdV'));
set(gca,'YDir','normal');
caxis ([-.03 .03]);
%title('dI/dV');
colormap(spdIdV,bluewhitered(256))
c=colorbar('southoutside');
min_dIdV = min(dIdV(dIdV<0)); % min(dIdV(dIdV<0)) is a workaround
if(isempty(min_dIdV))
    min_dIdV = min(dIdV(:));
end
c.Limits = [min_dIdV  max(dIdV(:))]; 
c.Label.String = 'dI/dV';
c.Label.FontSize = 12; c.FontSize = 12;
c.TickDirection='out';
c.TickLength = 0.05;

sp2eI = subplot(1,5,2);
schottky = abs(current);
Si = abs(current_noise);
mS = max( max(schottky(:)), max(Si(:)));
ih=imagesc(gateVv,biasVv,real(schottky'));
set(gca,'YDir','normal');
caxis([0, mS]);
%title('|2eI|');
colormap(sp2eI,flipud(gray))
c=colorbar('southoutside');
c.Label.String = '|2eI|';
c.Label.FontSize = 12; c.FontSize = 12;
c.TickDirection='out';
c.TickLength = 0.05;

spSi = subplot(1,5,3);

ih=imagesc(gateVv,biasVv,real(Si'));
set(gca,'YDir','normal');
caxis([0, mS]);
colormap(spSi,flipud(gray))
c=colorbar('southoutside');
c.Label.String = 'S_I';
c.Label.FontSize = 12; c.FontSize = 12;
c.TickDirection='out';
c.TickLength = 0.05;

spSep = subplot(1,5,4);
S_ep = Si - schottky;
ih=imagesc(gateVv,biasVv,real(S_ep'));
set(ih,'alphadata',~isnan(S_ep'))
set(gca,'YDir','normal');
caxis([min(S_ep(:)), max(S_ep(:))]);
%title('S_{EP}');
colormap(spSep,bluewhitered(256))
c=colorbar('southoutside');
c.Label.String = 'S_I^{EP}';
c.Label.FontSize = 12; c.FontSize = 12;
c.TickDirection='out';
c.TickLength = 0.05;

spF = subplot(1,5,5);
Fano_factor = Si./schottky; 
Fano_factor(schottky<1e-15)=nan;
Fano_factor = Fano_factor - 1;

ih=imagesc(gateVv,biasVv,real(Fano_factor'));
set(ih,'alphadata',~isnan(Fano_factor'))
set(gca,'YDir','normal');
% set(sp,'clim',[0 3]);
caxis([-1, max(8,max(Fano_factor(:)))]);
%title('F');
%colormap(spF,flipud(gray))
colormap(spF,bluewhitered(256))
c=colorbar('southoutside');
c.Limits = [-1 max(Fano_factor(:))];
c.Label.String = 'F-1';
c.Label.FontSize = 12; c.FontSize = 12;
c.TickDirection='out';
c.TickLength = 0.05;
%suptitle(config_dir); %in BioInformatics toolbox