% Two Bus Example, see: `Harmonic Optimal Power Flow with Transformer 
% Excitation` by F. Geth and T. Van Acker, pg. 7, § IV.A.

function mpc = two_bus_example
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone    Vmax	Vmin
mpc.bus = [
    1       4       0.000   0.000   0.00    0.00    1       1.00    0.00    155      1       1.10    0.90;
    2       1       100.000   0.000   0.00    0.00    1       1.00	0.00    35.4      1       1.10    0.90;
];

%% bus harmonic data 
%column_names%  nh_1    nh_2    nh_3    std
mpc.bus_harmonics = [
                1.000   0.000   0.000   'IEC61000-2-4:2002, Cl. 2'
                1.000   0.200   0.000   'IEC61000-2-4:2002, Cl. 2'
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     900.00    -900.00   1.05    100.0   1       999.00   -999.0;
];

%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    %1       2       0.004   0.008   0       10      10      10      1       0       1       -60     60;
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     gsh     r1      r2      vg      gnd1    gnd2    re1     xe1     re2         xe2         rate_a
mpc.xfmr = [
                1       2       0.1     0.0     0.000   0.000   'Yy0'  1       1       0.0     0.0     0.0         0.0         125; %0.000775 %0.135958171 % 0.0016864
                %1       2       0.13    0.07752 0.00211 0.00211 'Yy0'  1       1       0.0     0.0     0.0         0.0         125; 
                %3       4       0.0229  0.01420 0.00030 0.00030 'Yy0'   0       0       0.0     0.0     0.0         0.0         31.5;
                %5       6       0.0107  0.00270 0.00070 0.00070 'Dy11'  0       1       0.0     0.0     0.0000046   0.0000063   2.5;
                %7       8       0.0002  0.00100 0.00005 0.00005 'Dz0'   0       1       0.0     0.0     0.000025    0.000013    0.25; 
];