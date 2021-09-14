function [ resultado, velcom,velcomc ] = saltogoleiro(dat3d, lado , graficos)
%Analise do salto de goleiros em cobrancas de penaltis
%   Prof. Dr. Paulo Roberto Pereira Santiago - LaBioCoM/USP
%     14/03/2015
%
% Os pontos devem ser definidos na ordem:
% Primeira coluna com os quadros de imagem (frames) e as demais com as
% coordenadas x, y e z dos seguintes pontos anatomicos =
% 1-Cabeça distal, 2-Esterno, 3-Ombro Direito, 4-Cotovelo Direito,
% 5-Punho Direito, 6-Mao Direita, 7-Ombro Esquerdo, 8-Cotovelo Esquerdo,
% 10-Punho Esquerdo, 10-Mao Esquerda, 11-Crista Iliaca Direita,  
% 12-Joelho Direito, 13-Tornozelo Direito, 14-Calcanhar Direito, 
% 15-Halux Direito, 16-Crista Iliaca Esquerda, 17-Joelho Esquerdo, 
% 18-Tornozelo Esquerdo, 19-Calcanhar Esquerdo, 20-Halux Esquerdo, 
pkg load io
pkg load signal
pkg load image

if nargin == 2;graficos = 1;end
freq = 120;
dat = load(dat3d); 

close all
nl = size(dat,1);
temp = [(0:nl-1)/freq]'; % vetor tempo

dat = dat/1000;
dat = filtbutter(dat,4,120,'low');

[com_total] = COM3D(dat);


%%%% Retirada do pe do solo foot off

haluxD_Z0 = dat(:,46);
haluxE_Z0 = dat(:,61);

zeropadD = min(haluxD_Z0);
zeropadE = min(haluxE_Z0);

haluxD_Z = haluxD_Z0 - zeropadD;
haluxE_Z = haluxE_Z0 - zeropadE;

[row_tirapeE0] = find(haluxE_Z > 0.25);
[row_tirapeD0] = find(haluxD_Z > 0.25);

if size(row_tirapeE0,1) == 0
    [row_tirapeE0] = find(haluxE_Z > 0.15);
end
if  size(row_tirapeD0,1) == 0
    [row_tirapeD0] = find(haluxD_Z > 0.15);
end


if size(row_tirapeE0,1) == 0
    [row_tirapeE0] = find(haluxE_Z > 0.10);
end
if  size(row_tirapeD0,1) == 0
    [row_tirapeD0] = find(haluxD_Z > 0.10);
end

if size(row_tirapeE0,1) == 0
    [row_tirapeE0] = find(haluxE_Z > 0.05);
end
if  size(row_tirapeD0,1) == 0
    [row_tirapeD0] = find(haluxD_Z > 0.05);
end


if size(row_tirapeE0,1) == 0
    disp(' ')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('Nao tirou o pe esquerdo do solo!!!')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    row_tirapeE0 = 13;
end

if  size(row_tirapeD0,1) == 0
    disp(' ')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('Nao tirou o pe direito do solo!!!')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    
    row_tirapeD0 = 13;
end

row_tirapeE = row_tirapeE0(1) - 12;
row_tirapeD = row_tirapeD0(1) - 12;

TR_pe_tempD = temp(row_tirapeD);
    
TR_pe_tempE = temp(row_tirapeE);



%%%%% angulos do joelho Dir. Esq. 
crid = dat(:,32:34);
joed = dat(:,35:37);
tord = dat(:,38:40);
crie = dat(:,47:49);
joee = dat(:,50:52);
tore = dat(:,53:55);

VcoxD = crid - joed;
VcoxE = crie - joee;

VperD = tord - joed;
VperE = tore - joee;

ang_joe_D = NaN(nl,1);
ang_joe_E = NaN(nl,1);

for i = 1:nl
    
    VcoxDn = VcoxD(i,:) / norm(VcoxD(i,:));
    VcoxEn = VcoxE(i,:) / norm(VcoxE(i,:));

    VperDn = VperD(i,:) / norm(VperD(i,:));
    VperEn = VperE(i,:) / norm(VperE(i,:));

    ang_joe_D(i,:) = 180-(acos(dot(VcoxDn,VperDn)) * (180/pi));
    ang_joe_E(i,:) = 180-(acos(dot(VcoxEn,VperEn)) * (180/pi));
    
end

%%%%%%%%%%%%%%% Identificando o primeiro movimento de joelho %%%%%%%%%%%%%
[locs_ang_D] = startmov(ang_joe_D,round(freq/10),3); 
[locs_ang_E] = startmov(ang_joe_E,round(freq/10),3);


%%
%%%%%%%%%%%%% Velocidades %%%%%%%%%%%
% if lado == 'd'
%     com_total(:,2) = com_total(:,2) * -1;
% end


diffcom = diff(com_total);
normcom = NaN(size(diffcom,1),1);

for i = 1:size(diffcom,1)
    normcom(i,1) = norm(diffcom(i,:));
end

velcom = normcom / (1/freq);

velcom = [NaN;velcom];

vel_com_xyz = diffcom / (1/freq);


vel_com_xyz = [[NaN,NaN,NaN];vel_com_xyz];


[row_3D,col_3D] = find(velcom >= 0.9); % velocidade 3D

if lado == 'e'
[row,col] = find(vel_com_xyz(:,2) >= 0.9); % velocidade ml
end

if lado == 'd'
[row,col] = find(vel_com_xyz(:,2) <= -0.9); % velocidade ml
end


[maxvel,locmax] = max(velcom);
% [pks,locs] = findpeaks(velcom,'SORTSTR','descend');


diffvelcom_xyz = diff(vel_com_xyz);
diffvelcom = diff(velcom);

acelcom = diffvelcom_xyz ./ (1/freq);
acelcom_3D = diffvelcom ./ (1/freq);

acelcom = [[NaN,NaN,NaN];acelcom];
acelcom_3D = [NaN;acelcom_3D];

% [locs,datc] = startmov(acelcom,10,5)
if lado == 'e'
[row_acc,col_acc] = find(acelcom(:,2) >= 3); % aceleracao
end

if lado == 'd'
[row_acc,col_acc] = find(acelcom(:,2) <= -3); % aceleracao
end

if row(1) < row_acc(1)
    row_sel = row(1);
end


if row_acc(1) < row(1)
    row_sel = row_acc(1);
end

%%%%%%%% Angulo de saida
vet_saida = com_total(locmax,:) - com_total(row_sel,:);
    [azimuth,elevation,r] = cart2sph(vet_saida(:,2),vet_saida(:,1),vet_saida(:,3));
if lado == 'e'
    [ang_vet] = [rad2deg(azimuth),rad2deg(elevation),r];
end

if lado == 'd'&& azimuth < 0
    [ang_vet] = [180-(rad2deg(azimuth) + 360), rad2deg(elevation), r];
end

if lado == 'd'&& azimuth > 0
    [ang_vet] = [180-(rad2deg(azimuth)), rad2deg(elevation), r];
end


%%%%%%%%%%%% Resultaldos %%%%%%%%%%%%%%%%%5

TR = temp(row_sel); % tempo de resposta TR
tmax = temp(locmax); % tempo para atingir o pico de velocidade
velTR = velcom(row_sel); % velocidade em TR
tmaxTR = tmax - TR; % tempo para atigir o pico de velocidade ap�s TR
taxavel = maxvel - velTR; % tempo para atigir o pico de velocidade ap�s TR
inddesemp = taxavel/tmaxTR; % indicie de desempenho;
vel_com_TR_peD = velcom(row_tirapeD); % velocidade na hora de retirada do p� direito
vel_com_TR_peE = velcom(row_tirapeE); % velocidade na hora de retirada do p� esquerdo
velcomc = velcom(row_sel:locmax,:); % veliciade corte a partir do TR da vel do com
TR_joe_D = temp(locs_ang_D); 
TR_joe_E = temp(locs_ang_E); 


if graficos == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % % % % % % % % % % % % % % % % % % Graficos
    figure 
    plot3 (com_total(:,1),com_total(:,2),com_total(:,3),'.')
    hold on
    plot3 (com_total(end,1),com_total(end,2),com_total(end,3),'*r')
    plot3 (com_total(1,1),com_total(1,2),com_total(1,3),'*g')
    plot3 (com_total(row_sel,1),com_total(row_sel,2),com_total(row_sel,3),'*k','markersize',5)
    zlim([0 2.44])
    ylim([0 7.32])
    xlim([0 1.5])
    if lado == 'e'
        view([85,10])
    elseif lado == 'd'
        view([110,10])
    end
    daspect([1 1 1])
    box on
    xlabel('X (antero-posterior)')
    ylabel('Y (medio-lateral)')
    zlabel('Z (vertical)')
    grid on
    line([0 1],[0 0],[0 0],'color','red','linewidth',4)
    line([0 0],[0 1],[0 0],'color','green','linewidth',4)
    line([0 0],[0 0],[0 1],'color','blue','linewidth',4)
    rotate3d on
    
    figure
    subplot(3,1,1)
    plot(temp,com_total(:,1))
    hold on
    plot(temp(row_sel),com_total(row_sel,1),'>g')
    plot(temp(locmax),com_total(locmax,1),'or')
    title('Coordenadas 3D do COM')
    ylabel('X - ant/post')

    subplot(3,1,2)
    plot(temp,com_total(:,2))
    hold on
    plot(temp(row_sel),com_total(row_sel,2),'>g')
    plot(temp(locmax),com_total(locmax,2),'or')
    ylabel('Y - med/lat')

    subplot(3,1,3)
    plot(temp,com_total(:,3))
    hold on
    plot(temp(row_sel),com_total(row_sel,3),'>g')
    plot(temp(locmax),com_total(locmax,3),'or')
    ylabel('Z vertical')



    figure
    subplot(4,1,1)
    plot(temp,vel_com_xyz(:,1))
    hold on
    plot(temp(row_sel),vel_com_xyz(row_sel,1),'>g')
    plot(temp(locmax),vel_com_xyz(locmax,1),'or')
    title('Velocidades X, Y e Z do COM')
    ylabel('X - ant/post')

    subplot(4,1,2)
    plot(temp,vel_com_xyz(:,2))
    hold on
    plot(temp(row_sel),vel_com_xyz(row_sel,2),'>g')
    plot(temp(locmax),vel_com_xyz(locmax,2),'or')
    ylabel('Y - med/lat')

    subplot(4,1,3)
    plot(temp,vel_com_xyz(:,3))
    hold on
    plot(temp(row_sel),vel_com_xyz(row_sel,3),'>g')
    plot(temp(locmax),vel_com_xyz(locmax,3),'or')
    ylabel('Z - vertical')


    subplot(4,1,4)
    plot(temp,acelcom(:,2))
    hold on
    plot(temp(row_sel),acelcom(row_sel,2),'>g')
    plot(temp(locmax),acelcom(locmax,2),'or')
    title('Aceleracao medio-lateral do COM')
    ylabel('y aceleracao med-lat')
    xlabel('Tempo (s)')



    figure
    ax(1)=subplot(2,1,1);
    plot(temp,velcom,'k')
    title('Velocidade e Aceleracao do CoM')
    hold on
    h1 = plot(temp(row_sel(1)),velcom(row_sel(1)),'>g','markersize',6);
    h2 = plot(temp(locmax),maxvel,'or','markersize',6);
    h3 = plot(temp(row_tirapeD),velcom(row_tirapeD),'+r','markersize',6);
    h4 = plot(temp(row_tirapeE),velcom(row_tirapeE),'+b','markersize',6);

    if lado == 'd' 
        h101 = plot(temp(locs_ang_D),velcom(locs_ang_D),'sg','markersize',6);
    end
    if lado == 'e'
        h101 = plot(temp(locs_ang_E),velcom(locs_ang_E),'sg','markersize',6);
    end

    legend([h1,h101,h2,h3,h4],'TR CG','TR joelho','Vel. Max. CG','foot off Dir.','foot off Esq.')

    ylabel('velocidade [m/s]')

    hold off

    ax(2)=subplot(2,1,2);
    plot(temp,acelcom_3D,'k')
    hold on
    plot(temp(row_sel),acelcom_3D(row_sel),'>g','markersize',6)
    plot(temp(locmax),acelcom_3D(locmax),'or','markersize',6);
    plot(temp(row_tirapeD),acelcom_3D(row_tirapeD),'+r','markersize',6);
    plot(temp(row_tirapeE),acelcom_3D(row_tirapeE),'+b','markersize',6);
    if lado == 'd' 
        plot(temp(locs_ang_D),acelcom_3D(locs_ang_D),'sg','markersize',6);
    end
    if lado == 'e'
        plot(temp(locs_ang_E),acelcom_3D(locs_ang_E),'sg','markersize',6);
    end

    xlabel('tempo [s]')
    ylabel('aceleracao [m/s^2]')

    linkaxes([ax(1),ax(2)],'x');



    figure
    plot(temp,haluxE_Z,'b',temp,haluxD_Z,'r')
    legend('esquerdo','direito')
    hold on
    plot(temp(row_tirapeE),haluxE_Z(row_tirapeE),'+b',temp(row_tirapeD),haluxD_Z(row_tirapeD),'+r')
    title('Coordenada vertical do halux')
    ylabel('metros')
    xlabel('segundos')


    figure
    plot(temp,ang_joe_E,'b',temp,ang_joe_D,'r')
    legend('esquerdo','direito')
    hold on
    h7 = plot(temp(row_sel(1)),ang_joe_E(row_sel(1)),'>g','markersize',6);
    h8 = plot(temp(row_sel(1)),ang_joe_D(row_sel(1)),'>g','markersize',6);
    h5 = plot(temp(row_tirapeE),ang_joe_E(row_tirapeE),'+b');
    h6 = plot(temp(row_tirapeD),ang_joe_D(row_tirapeD),'+r');
    h9 = plot(temp(locmax),ang_joe_E(locmax),'ok','markersize',6);
    h10 = plot(temp(locmax),ang_joe_D(locmax),'ok','markersize',6);
    h11 = plot(temp(locs_ang_D),ang_joe_D(locs_ang_D),'sk');
    h12 = plot(temp(locs_ang_E),ang_joe_E(locs_ang_E),'sk');

    legend([h8,h5,h6,h9,h11],'TR','foot off E','foot off D','Vel. max.','TR joelho')
    title('ANGULO DO JOELHO')
    ylabel('graus')
    xlabel('segundos')

end








disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('   ')
disp('###########################################################################################################################')
disp('                    |      Coluna 1      |       Coluna 2      |      Coluna 3        |   Coluna 4      |    Coluna 5 ')
disp('--------------------|--------------------|---------------------|----------------------|-----------------|------------------')
disp('lin 1 (tempos)      |         TR         |   temp. max. vel.   | temp. TR a max. vel. |   TR_pe dir.    |    TR_pe esq.')
disp('lin 2 (velocidades) |     vel. em TR     |      max. vel.      |    ind. desempenho   | vel. TR_pe dir. | Vel. TR_pe esq.')
disp('lin 3 Angulos       | ang. saida frontal | ang. saida vertical |    dist  max. vel.   | TR joelho Dir.  |  TR joelho Esq.')
disp('--------------------|------------------------------------------------------------------------------------------------------')
disp('###########################################################################################################################')
disp('   ')

resultado = [TR , tmax , tmaxTR , TR_pe_tempD , TR_pe_tempE ;...
             velTR , maxvel , inddesemp , vel_com_TR_peD , vel_com_TR_peE ;...
             ang_vet(1) , ang_vet(2) , ang_vet(3) , TR_joe_D , TR_joe_E]
         
         disp('#####################################################################################################')

         if TR < 0.1
             disp('GOLEIRO ANTECIPOU!!!!!!!')

         end
end






function [com_total] = COM3D(dat)
% CALCULO DO COM EM IMAGEM (TRIDIMENSIONAL)
% Rotina criada por Reinaldo, Natalia e Priscila.
% Primeira versao 15/06/13
%
% Disciplina: Metodologia Analise Cinem�tica (USP)
%
% [CG_total,CG_seg,] = com3d(X)
%
% Coordenadas tridimensionais (3D) 20 marcadores
% Os pontos devem ser definidos na ordem: 
% 1-Cabeca distal, 2-Externo, 3-Ombro Direito, 4-Cotovelo Direito,
% 5-Punho Direito, 6-Mao Direita, 7-Ombro Esquerdo, 8-Cotovelo Esquerdo,
% 10-Punho Esquerdo, 10-Mao Esquerda, 11-Crista Ilaaca Direita,  
% 12-Joelho Direito, 13-Tornozelo Direito, 14-Calcanhar Direito, 
% 15-Halux Direito, 16-Crista Ilaaca Esquerda, 17-Joelho Esquerdo, 
% 18-Tornozelo Esquerdo, 19-Calcanhar Esquerdo, 20-Halux Esquerdo, 
% 
%
% #########################################################################
% Tabela de massas relativas e localizacao dos centros de massa de cada
% segmento (Plagenhoef et al., 1983)
% #########################################################################
% Segmento     Massas relativas    localizacao COM      Segmento
%                                                    
% Cabeca+P.       0.0823               0.55         vertex - 7 cervical
% Tronco          0.46,01              0.4175       esterno - media crisIliacas
% Braco           0.03075              0.447        c. glenohumeral - art.cotovelo
% Antebraco       0.0172               0.432        art.cotovelo -    art.pulso
% Mao             0.00575              0.468        art.pulso -  art. 3 falange
% Coxa            0.11125              0.4305       crisIliac - cond. femural
% Perna           0.05050              0.4265       cond. femural - maleolo
% Pe              0.0138               0.5          halux - calcanio

% dat = filtbutter(dat,4,120,'low');
% Cabeca

cab = dat(:,2:4);
ext = dat(:,5:7);

com_cab = ext + 0.55 * (cab - ext);


% hold on
% segcab = [ore_d;ore_e];
% plot(segcab(:,1),segcab(:,2),segcab(:,3),'linewidth',3,'color','y')
% plot(ore_d(1),ore_e(2),'.','MarkerSize',16);
% plot(ore_d(1),ore_e(2),'.','MarkerSize',16);
% plot(com_cab(1),com_cab(2),'r.','MarkerSize',16);

% Tronco

tron_d = (dat(:,32:34) + dat(:,47:49))/2;
com_tr = ext + 0.4175 * (tron_d - ext);


% Bracos

ombd = dat(:,8:10);
cotd = dat(:,11:13);
ombe = dat(:,20:22);
cote = dat(:,23:25);

com_brd = ombd + 0.447 * (cotd - ombd);
com_bre = ombe + 0.447 * (cote - ombe);

% AnteBracos

pund = dat(:,14:16);
pune = dat(:,26:28);

com_abd = cotd + 0.432 * (pund - cotd);
com_abe = cote + 0.432 * (pune - cote);

% Maos

maod = dat(:,17:19);
maoe = dat(:,29:31);

com_mad = pund + 0.468 * (maod - pund);
com_mae = pune + 0.468 * (maoe - pune);

% Coxas 

crid = dat(:,32:34);
joed = dat(:,35:37);
crie = dat(:,47:49);
joee = dat(:,50:52);

com_cxd = crid + 0.4305 * (joed - crid);
com_cxe = crie + 0.4305 * (joee - crie);

% Pernas 

tord = dat(:,38:40);
tore = dat(:,53:55);

com_pnd = joed + 0.4265 * (tord - joed);
com_pne = joee + 0.4265 * (tore - joee);

% Pes

hald = dat(:,44:46);
cald = dat(:,41:43);
hale = dat(:,59:61);
cale = dat(:,56:58);

com_ped = hald + 0.500 * (cald - hald);
com_pee = hale + 0.500 * (cale - hale);


%Calculo do centro de massa total


com_total = ((0.0823 * com_cab) + (0.4601 * com_tr) + (0.03075 * com_brd) + (0.03075  * com_bre) + ...
            (0.0172 * com_abd) + (0.0172 * com_abe) + (0.00575 * com_mad) + (0.00575 * com_mae) + ...
            (0.11125 * com_cxd) + (0.11125 * com_cxe) + (0.0505 * com_pnd) + (0.0505 * com_pne) + ...
            (0.0138 * com_ped) + (0.0138 * com_pee)) / 1;
       
% com_seg = [com_cab , com_tr , com_brd  , com_bre , com_abd , com_abe , com_mad ,...
%          com_mae , com_cxd , com_cxe , com_pnd , com_pne , com_ped , com_pee]';
       
end



function [datf] = filtbutter(dat,fc,freq,ftype)
if nargin == 2; 
    freq = 100;
    ftype = 'low'; 
end
if nargin == 3; 
    ftype = 'low'; 
end

n=4; %ordem do filtro

wn=fc/(freq/2);   %frequencia de corte de

[b,a] = butter(n,wn,ftype); %definindo os parametros para o filtro de Butterworth

[nlin,ncol] = size(dat);
 
datf = NaN(nlin,ncol);
for i = 1:ncol;
datf(:,i) = filtfilt(b,a,dat(:,i));
end

end


function [inicio,datc] = startmov(dat,X,Y)
% selsplat(dat,X,Y)
%
% Rotina criada por Paulo R. P. Santiago (Preto).
% Primeira versao 24/11/06, ultima atualizacao ?????.
% Seleciona o sinal de inicio fim da forca vertical na fase durante
% o apoio na plataforma de forca (fz)
% dat = sinal da plataforma 
% 
% Y = limiar para o inicio do corte
%
% X = numero de pontos iniciais para definir o limiar de corte 

freq1 = 120;
if nargin==1, X = freq1; Y = 2;end
if nargin==2, Y = 2;end

nlin = size(dat,1);

% freq = X;

ddat = diff(dat);
% ddat = dat;

% pcini = Y * (max(ddat(1:freq)));
med_ddat = mean(ddat(1:X));
std_ddat = std(ddat(1:X));

pcini_u = med_ddat + (Y * std_ddat);
pcini_d = med_ddat - (Y * std_ddat);

pcini1_u = find(ddat > pcini_u);
pcini1_d = find(ddat < pcini_d);

if size(pcini1_u,1) == 0
    pcini1_u = 1;
end

if size(pcini1_d,1) == 0
    pcini1_d = 1;
end

if  pcini1_u(1) < pcini1_d(1)
    inicio = pcini1_u(1);
end

if  pcini1_u(1) > pcini1_d(1)
    inicio = pcini1_d(1);
end


datc = dat(inicio:end,:);


end
