% CONFIG: Mixing pulse at variable time-delay since BEC collision
% DKS
% 2019-02-04

% general
exp_param=expparams();
vz=exp_param.vz;

% data import
configs.path.base='C:\Users\HE BEC\data\loop_tdelayscan_mf1';
configs.path.data=configs.path.base;
configs.load.path=fullfile(configs.path.data,'d');

configs.flag.param_scan=1;      % 0 for no scan; 1 for param scan
configs.path.paramlog=fullfile(configs.path.data,'LOG_parameters.txt');

configs.load.id=[];
configs.load.mincount=[];
configs.load.maxcount=[];     

configs.load.window{1}=[0.37,0.44];         % T [s]
configs.load.window{2}=[-35e-3,35e-3];      % X [m]
configs.load.window{3}=[-35e-3,35e-3];      % Y [m]

% halo
configs.mf(1).mf=1;
configs.mf(1).window={vz*[0.3822,0.402],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(1).p_bec0=[vz*0.3866,  -2.8e-3,    3.2e-3;
                   vz*0.3993,  -3.5e-3,    3.2e-3];
configs.mf(1).r_bec0=5e-3;

configs.mf(2).mf=0;
configs.mf(2).window={vz*[0.402,0.422],35e-3*[-1,1],35e-3*[-1,1]};
configs.mf(2).p_bec0=[vz*0.4051,  0.9e-3,    3.2e-3;
                    vz*0.4181,  1.1e-3,    3.2e-3];
configs.mf(2).r_bec0=5e-3;
                

% filter
% 1st stage
configs.filt.r_ball=12e-3;
configs.filt.r_crop=[0.6,1.2];
configs.filt.z_cap=0.85;

% post-distortion
configs.filt2.r_crop=[0.9,1.1];
configs.filt2.z_cap=0.75;


% analysis
configs.roi = 'all';       % 'halo', 'all'