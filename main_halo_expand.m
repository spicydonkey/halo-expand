% Mixing pulse at variable time-delay since BEC collision
%   effects: halo formation, free-fall, beam profile
% DKS
% 2019-02-04

%% configs
config_name='C:\Users\HE BEC\Documents\MATLAB\halo-expand\config_v1.m';

% VIS
config_fig.units='centimeters';
config_fig.pos_full=[0,0,17.2,6];       % full page width
config_fig.pos_2col=[0,0,8.6,3.2];      % 2-column
config_fig.rend='painters';
config_fig.ax_fontsize=9;
config_fig.ax_lwid=0.5;
config_fig.mark_siz=4;
config_fig.line_wid=0.5;


%% load config
run(config_name);


%% Get experimental params
% load wfmgen log
[params,id_in_param,param_id,Ipar]=wfmgen_log_parser(configs.path.paramlog);
nparam=size(params,1);      % number of unique param-set

% manually fix bug in LOG: file IDs skipped 1
id_in_param=cellfun(@(x) x+1,id_in_param,'UniformOutput',false);

t_exp=params(:,3);      % get expansion times (s)
nt=numel(t_exp);

%% load txy
[txy,fout]=load_txy(configs.load.path,configs.load.id,configs.load.window,configs.load.mincount,configs.load.maxcount,[],1,2,0);

zxy=txy2zxy(txy,vz);
shot_id=fout.id_ok;    
% ensure column vector
if size(shot_id,2)>1
    shot_id=shot_id';
end

% build boolean array selector for scanned param-set
b_paramset=cellfun(@(idx) ismember(shot_id,idx),...
    id_in_param,'UniformOutput',false);
b_paramset=horzcat(b_paramset{:});


% DEBUG
h_zxy_raw=figure;
plot_zxy(zxy,1e5);
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%% distinguish mJ and capture halo
n_shot=size(zxy,1);
n_mf=numel(configs.mf);

% preallocate
p_bec=cell(1,n_mf);     % marker BEC center positions
p_bec0=cell(1,n_mf);    % marker BEC in halo centered ZXY coord sys
p_halo=cell(1,n_mf);   % approx halo center (mid-point of marker BECs)
for ii=1:n_mf
    p_bec{ii}=NaN(n_shot,3,2);
    p_halo{ii}=NaN(n_shot,3);
end
zxy0=cell(n_shot,n_mf); % halo centerised zxy


for ii=1:n_shot
    tzxy=zxy{ii};
    for jj=1:n_mf
        tzxy_mf=boxcull(tzxy,configs.mf(jj).window);
        tp_bec0=configs.mf(jj).p_bec0;
        tr_bec0=configs.mf(jj).r_bec0;
        
        for kk=1:size(tp_bec0,1)
            tp_bec=capture_bec(tzxy_mf,tp_bec0(kk,:),tr_bec0,0);
            p_bec{jj}(ii,:,kk)=tp_bec;
        end
        tp_halo=mean(p_bec{jj}(ii,:,:),3);
        p_halo{jj}(ii,:)=tp_halo;
        zxy0{ii,jj}=tzxy_mf-tp_halo;
        p_bec0{jj}(ii,:,:)=p_bec{jj}(ii,:,:)-tp_halo;
    end
end

% % DEBUG
% figure(h_zxy_raw);
% hold on;
% plot_zxy(p_halo,[],50,'mk');

h_zxy0=figure();
plot_zxy(zxy0,3e5);
hold on;
for ii=1:2
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%% filter data
zxy0_filt=zxy0;     % initialise filtered data

%%% BEC + thermal
r_thermal=configs.filt.r_ball;
for ii=1:n_shot
    for jj=1:n_mf
        for kk=1:2      % pair of marker BECs
            tp_bec0=p_bec0{jj}(ii,:,kk);
            % get atoms outside ball centered at BEC
            zxy0_filt{ii,jj}=cropBall(zxy0_filt{ii,jj},r_thermal,tp_bec0,false);
        end
    end
end

% DEBUG
h_zxy0_filt_1=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%%% radial 
% here we rely on "average radius" of halo determined by the mean marker
% BEC locations

% estimate halo radius from marker BECs
r_crop=configs.filt.r_crop;
r_halo=NaN(n_shot,n_mf);
for ii=1:n_mf
    r_halo(:,ii)=0.5*vnorm(p_bec0{ii}(:,:,2)-p_bec0{ii}(:,:,1));
end
r_halo_avg=mean(r_halo,1);

% filter
for ii=1:n_mf
    tr_lim=r_crop*r_halo_avg(ii);       % absolute norm limits
	zxy0_filt(:,ii)=cfilter_norm(zxy0_filt(:,ii),tr_lim(1),tr_lim(2));
end

% DEBUG
h_zxy0_filt_2=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);


%%% polar caps
% build box window for halo between caps
z_cap=configs.filt.z_cap;
window_z_filt=cell(1,n_mf);
for ii=1:n_mf
    window_z_filt{ii}={z_cap*r_halo_avg(ii)*[-1,1],[],[]};
end

% filter
for ii=1:n_mf
    tbox_lim=window_z_filt{ii};
    zxy0_filt(:,ii)=cellfun(@(x) boxcull(x,tbox_lim),zxy0_filt(:,ii),'UniformOutput',false);
end

% DEBUG
h_zxy0_filt_3=figure();
plot_zxy(zxy0_filt,3e5);
hold on;
for ii=1:n_mf
    for jj=1:2
        plot_zxy({p_bec0{ii}(:,:,jj)},[],100,'k');
    end
end
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view([0,0]);



%% k-space and distortion cancellation
k_halo=zxy0_filt;       % initialise atoms in k-space (i.e. atoms lie on unit-sphere)
v_ellip=cell(n_mf,1);

% ellipsoid fit to sphere
for ii=1:n_mf
    [k_halo(:,ii),v_ellip{ii}]=map2usph(k_halo(:,ii));
end
v_ellip=[v_ellip{:}];   % form as struct array
    

% VIS
scatter_halo(k_halo);


%% transform raw spatial clouds
k_all=cell(size(zxy0));
for ii = 1:n_mf
    tV = v_ellip(ii);       
    k_all(:,ii) = cellfun(@(x) ellip2usph(x,tV.cent,tV.rad,tV.vaxis),zxy0(:,ii),...
        'UniformOutput',false);
end

% VIS
scatter_halo(k_all);

%%% filter radially
k_all_filt = cfilter_norm(k_all,configs.filt2.r_crop(1),configs.filt2.r_crop(2));

% VIS
scatter_halo(k_all_filt);


%% filter post-processed data
%%% radial
r_crop2=configs.filt2.r_crop;
k_halo_filt=cfilter_norm(k_halo,r_crop2(1),r_crop2(2));

%%% z-cap
% build box window for halo between caps
z_cap2=configs.filt2.z_cap;
window_z_filt2=cell(1,n_mf);
for ii=1:n_mf
    window_z_filt2{ii}={z_cap2*[-1,1],[],[]};
end

% filter
for ii=1:n_mf
    tbox_lim=window_z_filt2{ii};
    k_halo_filt(:,ii)=cellfun(@(x) boxcull(x,tbox_lim),k_halo_filt(:,ii),'UniformOutput',false);
end


% DEBUG
scatter_halo(k_halo_filt);


%% categorise data by exp-params
k_par=cell(1,nparam);
if nparam>1
    for ii=1:nparam
        k_par{ii}=k_halo_filt(b_paramset(:,ii),:);      % get all halo data
        %from this param-set and store
    end
else
    k_par{1}=k_halo_filt;
end

n_shot_par=shotSize(k_par);

% DEBUG
figure;
for ii=1:nparam
    subplot(1,nparam,ii);
    plot_zxy(k_par{ii});
    
    title(num2str(1e3*t_exp(ii),2));
    axis equal;
    xlabel('$k_x$');
    ylabel('$k_y$');
    zlabel('$k_z$');
    view([0,0]);
end

%%% END OF PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PRE-analysis
% SAVE ORIGINAL 
if ~exist('k_par_orig','var')
    k_par_orig=k_par;
end

%%% 3D rotated orientation
% v_euler=[-pi/2,pi*3/4,-pi/2];             % proper euler angle
% invR=inv(euler2rotm2(v_euler));     % rotation matrix: xyz --> XYZ

v_euler=[0,0,0];                    % original coord
invR=inv(euler2rotm2(v_euler));     % rotation matrix: xyz --> XYZ

% transform to RH coords
% k_par_rh=cellfun(@(C) cellfun(@(xl) tzxy2RHtzxy(xl),C,'UniformOutput',false),...
%     k_par_orig,'UniformOutput',false);    % z along g
k_par_rh=cellfun(@(C) cellfun(@(xl) tzxy2RHtzxy2(xl),C,'UniformOutput',false),...
    k_par_orig,'UniformOutput',false);      % EXP-coord sys (z against g)

% transform to xyz coords
k_par_xyz=cellfun(@(C) cellfun(@(zxy) zxy2xyz(zxy),C,'UniformOutput',false),...
    k_par_rh,'UniformOutput',false);

% transform to XYZ coords (rotated)
k_par_XYZ=cellfun(@(C) cellfun(@(xyz) (invR*xyz')',C,'UniformOutput',false),...
    k_par_xyz,'UniformOutput',false);

% transform back to ZXY coords
k_par=cellfun(@(C) cellfun(@(xyz) xyz2zxy(xyz),C,'UniformOutput',false),...
    k_par_XYZ,'UniformOutput',false);


%% ANALYSIS
%%% configs
% momentum zones
configs.zone.lim_az=pi*[-1,1];
configs.zone.lim_el=pi/2*[-1,1];        % pi/4

configs.zone.n_az=120;
configs.zone.n_el=60;
configs.zone.alpha=pi/20;       % pi/10 half-cone angle

% % DEBUG
% configs.zone.n_az=20;
% configs.zone.n_el=10;
% configs.zone.alpha=pi/10;

configs.zone.az=linspace(configs.zone.lim_az(1),configs.zone.lim_az(2),configs.zone.n_az+1);
configs.zone.az=configs.zone.az(1:end-1);
configs.zone.el=linspace(configs.zone.lim_el(1),configs.zone.lim_el(2),configs.zone.n_el);

azel=configs.zone;

[az_grid,el_grid]=ndgrid(configs.zone.az,configs.zone.el);    % az-el grid
n_zone=numel(az_grid);
azel_dims=size(az_grid);

%% momentum-spin resolved atom numbers
Nm_k=cell(nt,1);     % #atoms in k,mj-zone categorise by exp param

% KS counting
for ii=1:nt
    tk=k_par{ii};        % exp data for this param
    
    % get num in zone/shot
    tN=arrayfun(@(th,ph) cellfun(@(x) size(inCone(x,th,ph,azel.alpha),1),tk),...
        az_grid,el_grid,'UniformOutput',false);     
    
    % format into multi-dim array: AZ X EL X SHOT X MJ
    ttN=NaN([azel_dims,n_shot_par(ii),n_mf]);
    for jj=1:n_zone
        [iaz,iel]=ind2sub(azel_dims,jj);
        ttN(iaz,iel,:,:)=tN{iaz,iel};
    end
    Nm_k{ii}=ttN;
end

M_k=cellfun(@(x) -diff(x,1,4),Nm_k,'UniformOutput',false);  % polarisation/(M)agnetisation
N_k=cellfun(@(x) sum(x,4),Nm_k,'UniformOutput',false);      % total number (N)

% normalise
m_k=cellfun(@(x,n) x./n,M_k,N_k,'UniformOutput',false);

m_k_avg=cellfun(@(x) mean(x,3,'omitnan'),m_k,'UniformOutput',false);
m_k_avg=cat(3,m_k_avg{:});
m_k_std=cellfun(@(x) std(x,0,3,'omitnan'),m_k,'UniformOutput',false);
m_k_std=cat(3,m_k_std{:});
m_k_se=m_k_std./sqrt(reshape(n_shot_par,[1,1,nt]));


%% CHECK: experiment stability
%%% config
% col_texp=gray(nt+1);
col_texp=parula(nt+1);

% regions to display
nkdisp=4;
disp_iazel=round(linspace(1,configs.zone.n_az,nkdisp+1))';      % uniform around equator [0,pi]
disp_iazel=disp_iazel(1:end-1);
disp_iazel(:,2)=round(configs.zone.n_el/2);

%%% vis
figure('Name','exp_stability');
for jj=1:nkdisp
    subplot(1,nkdisp,jj);
    hold on;
    for ii=1:nt
        plot(squeeze(m_k{ii}(disp_iazel(jj,1),disp_iazel(jj,2),:)),'o-',...
            'Color',col_texp(ii,:));
        
        xlabel('shot \#');
        ylabel('polarisation')
        ylim([-1,1]);
    end
end

%% VIS: Momentum distribution of polarisation
figname='pol_kdist_vs_texp';
h=figure('Name',figname,'Units','centimeters','Position',[0,0,13,38],'Renderer','opengl');

for ii=1:nt
% for ii=1:5
    subplot(nt,1,ii);
    
    tmap=plotFlatMapWrappedRad(az_grid,el_grid,m_k_avg(:,:,ii),'rect','texturemap');
    
    ax=gca;
    axis tight;
    xlim([-180,180]);
    xticks(-180:90:180);
    yticks(-45:45:45);
    xlabel('$\theta$ (deg)');
    ylabel('$\phi$ (deg)');
    axis off;
    
    caxis([-1,1]);
%     colormap('parula');       % discrete bands at the slope
    colormap('magma');
    
    temp_str=sprintf('%0.2g %s',1e3*t_exp(ii),'ms');
    text(-210,0,temp_str,'Rotation',-90,'HorizontalAlignment','center',...
        'VerticalAlignment','middle','FontSize',12)
    
    if ii==1
        cbar=colorbar('north');
        cbar.TickLabelInterpreter='latex';
        cbar.Label.Interpreter='latex';
        cbar.Label.String='polarisation';
        cbar.Label.FontSize=12;
        cbar.FontSize=12;
        
        %%% manual ticks
        xlim0=get(ax,'XLim');
        ylim0=get(ax,'YLim');
        text(mean(xlim0),max(ylim0)+30,'$\theta$ (deg)',...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
        text(max(xlim0)+30,mean(ylim0),'$\phi$ (deg)',...
                'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10,...
                'Rotation',-90);
        for xx=xticks
            text(xx,max(ylim0)+10,num2str(xx),...
                'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10); 
        end
        for yy=yticks
            text(max(xlim0)+10,yy,num2str(yy),...
                'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10);
        end
    end
end
% locate colorbar to bottom
%   NOTE: current axes is to the LAST (bottom) subplot)
pos_bottom=ax.Position;
% define colorbar dims relative to subplot axes
cbar_width=pos_bottom(3);           
cbar_height=pos_bottom(4)/6;         
cbar_x=pos_bottom(1);
cbar_y=pos_bottom(2)-2.5*cbar_height;
pos_cbar=[cbar_x, cbar_y, cbar_width, cbar_height];

cbar.Position=pos_cbar;
cbar.AxisLocation='out';        % label outward

% annotation
AxesH = axes('Parent', h, ...
  'Units', 'normalized', ...
  'Position', [0, 0, 1, 1], ...
  'Visible', 'off', ...
  'XLim', [0, 1], ...
  'YLim', [0, 1], ...
  'NextPlot', 'add');
text('parent',AxesH,'Units','normalized','Position',[0,0.5],...
    'Rotation',-90,'String','time since collision',...
    'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom');

% to print:
% print_vecrast(h,'exportType','svg');

%% VIS: Polarisation vs. t-expansion
col_nkdisp={'b','r','k','g'};
colf_nkdisp={'w','w','k','g'};
mrk_nkdisp={'o','^','s','d'};


figname='pol_vs_texp';
h=figure('Name',figname,'Units',config_fig.units,'Position',[0,0,8.6,6.5],...
    'Renderer',config_fig.rend);
hold on;

pleg=NaN(nkdisp,1);
for ii=1:nkdisp
    tp=ploterr(1e3*t_exp,squeeze(m_k_avg(disp_iazel(ii,1),disp_iazel(ii,2),:)),...
        [],squeeze(m_k_std(disp_iazel(ii,1),disp_iazel(ii,2),:)),'or-');
    set(tp(1),'Color',col_nkdisp{ii},'MarkerFaceColor',colf_nkdisp{ii},'Marker',mrk_nkdisp{ii},...
        'MarkerSize',config_fig.mark_siz,'LineWidth',config_fig.line_wid);
    set(tp(2),'Color',col_nkdisp{ii},'LineWidth',config_fig.line_wid);
    
    taz=configs.zone.az(disp_iazel(ii,1));
    tel=configs.zone.el(disp_iazel(ii,2));
    tstr_azel=sprintf('%3.2g, %3.2g',rad2deg(taz),rad2deg(tel));
    set(tp(1),'DisplayName',tstr_azel);
    
    pleg(ii)=tp(1);
end

ax=gca;
box on;
set(ax,'Layer','Top');
set(ax,'FontSize',config_fig.ax_fontsize);
set(ax,'LineWidth',config_fig.ax_lwid);

xlabel('delay (ms)');
ylabel('polarisation');

xlim([0,9]);
ylim(1.1*[-1,1]);

lgd=legend(pleg,'Location','SouthEast');
legend boxoff
lgd.Title.String='$\theta,\phi$ (deg)';
uistack(lgd,'bottom')