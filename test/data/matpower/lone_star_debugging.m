% Industrial Network, see: `Harmonic Optimal Power Flow with Transformer
% Excitation` by F. Geth and T. Van Acker, pg. 8, § IV.B.

function mpc = industrial_network
mpc.version = '2';
mpc.baseMVA =  100.00;
mpc.principle = 'maximin'

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    0       4       0.000   0.000   0.00    0.00    1       1.00    0.00    155.0   1       1.00    1.00;
    1       1       0.000   0.000   0.00    0.00    1       1.00    0.00    155.0   1       1.06    0.94;
    2       1       50.00   15.88   0.00    0.00    1       1.00	0.00    35.4    1       1.06    0.94;
    3       1       0.000   0.000   0.00    0.00    1       1.00    0.00    35.4    1       1.06    0.94;
    4       1       12.60   4.000   0.00    0.00    1       1.00	0.00    10.0    1       1.10    0.90;
    5       1       0.000   0.000   0.00    0.00    1       1.00    0.00    10.0    1       1.10    0.90;
    6       1       1.420   0.000   0.00    0.00    1       1.00	0.00    0.69    1       1.10    0.90;
    7       1       0.000   0.000   0.00    0.00    1       1.00    0.00    0.69    1       1.10    0.90;
    8       1       0.050   0.040   0.00    0.00    1       1.00	0.00    0.40    1       1.10    0.90;
];

%% bus harmonic data 
%column_names%  ref_angle   std 
mpc.bus_harmonics = [
                0.0         'Clean Bus';                % 0
                0.0         'IEC61000-3-6:2008';        % 1
                -30.0       'IEC61000-2-4:2002, Cl. 2'; % 2
                -30.0       'IEC61000-2-4:2002, Cl. 2'; % 3
                -30.0       'IEC61000-2-4:2002, Cl. 2'; % 4
                -30.0       'IEC61000-2-4:2002, Cl. 2'; % 5
                -60.0       'IEC61000-2-4:2002, Cl. 2'; % 6
                -60.0       'IEC61000-2-4:2002, Cl. 2'; % 7
                -60.0       'IEC61000-2-4:2002, Cl. 2'; % 8
]

%% branch data
%   f_bus	t_bus	r	        x	        b	        rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    0       1       0.0024      0.0476      0.00000     1000    1000    1000    1       0       1       -60     60; % Ssc=2.1e9, XRr=20     
    2       3       0.0127      0.0134      0.00160     125     1000    1000    1       0       1       -60     60; % NA2XS2Y 3x1x400/35, l=1.5 km, r=0.106Ω/km, C=0.271 uF/km, L=0.357 mH/km 
    4       5       0.0332      0.0212      0.00003     31.5    100     100     1       0       1       -60     60; % NA2XS2Y 3x1x240/25, l=0.2 km, r=0.166Ω/km, C=0.456 uF/km, L=0.338 mH/km
    6       7       1.3615      1.1317      0.00000     2.5     100     100     1       0       1       -60     60; % NYY 3x240/120, l=0.07 km, r=0.0926Ω/km, L=0.245mH/km
];

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    0       0.0     0.0     250.00  -250.00 1.00    100.0   1       250.00  0;
];
%column_names%  inf     gsc     bsc xr_ratio
mpc.gen_imp = [
                1       0.0     0.0 30.0;
]

%% transformer data
%column_names%  f_bus   t_bus   xsc         gsh         r1          r2          vg      gnd1    gnd2    re1     xe1     re2         xe2         rate_a
mpc.xfmr = [
                1       2       0.136       0.000775    0.0017      0.0017      'Yd11'  0       0       0.0     0.0     0.0         0.0         125; 
                3       4       0.400       0.000142    0.0050      0.0050      'Yy0'   0       0       0.0     0.0     0.0         0.0         31.5;
                5       6       2.400       0.000027    0.1600      0.1600      'Dy11'  0       1       0.0     0.0     0.0000046   0.0000063   2.5;
                7       8       10.00       0.000010    1.2800      1.2800      'Dz0'   0       1       0.0     0.0     0.000025    0.000013    0.25; 
];