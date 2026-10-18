%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   5-node star test case (recombination isolation experiment).          %
%                                                                       %
%     node 2   node 3   node 4   node 5     (loads / harmonic units)     %
%        \       |        |       /                                      %
%         \______|__node 1__|_____/          (gen, type 3, weak grid)    %
%                                                                       %
%   Four IDENTICAL branches 1-2,1-3,1-4,1-5.  The generator bus reaches  %
%   harmonic ground through a LARGE short-circuit impedance, so its       %
%   harmonic voltage = Z_gen*(I2+I3+I4+I5) aggregates all four units.     %
%   No line charging ⇒ no resonance.                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = case5_shhc
mpc.version = '2';
mpc.baseMVA = 100.0;
mpc.principle = 'absolute equality'

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 0.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
	2	 1	 10.0	 5.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
	3	 1	 10.0	 5.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
	4	 1	 10.0	 5.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
	5	 1	 10.0	 5.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
];

%% bus harmonic data
%column_names%  ref_angle   standard
mpc.bus_harmonics = [
                0.0         'IEC61000-3-6:2008';       % 1  (gen / ref)
                0.0         'IEC61000-3-6:2008';       % 2  (harmonic unit)
                0.0         'IEC61000-3-6:2008';       % 3  (harmonic unit)
                0.0         'IEC61000-3-6:2008';       % 4  (harmonic unit)
                0.0         'IEC61000-3-6:2008';       % 5  (harmonic unit)
]

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
	1	 40.0	 20.0	 100.0	 -100.0	 1.0	 100.0	 1	 200.0	 0.0; % grid
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	 0.0	 0.0	 0.0; % grid
];

%% generator short-circuit impedance data (UNCHANGED: weak grid)
%column_names%  rsc     xsc
mpc.gen_sc = [
	0.10	2.00;   % large short-circuit impedance (gen bus aggregates)
]

%% branch data (four identical lines; b = 0)
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 2	 0.0002	 0.002	 0.0	 100	 100	 100	 0.0	 0.0	 1	 -360	 360;
	1	 3	 0.0002	 0.002	 0.0	 100	 100	 100	 0.0	 0.0	 1	 -360	 360;
	1	 4	 0.0002	 0.002	 0.0	 100	 100	 100	 0.0	 0.0	 1	 -360	 360;
	1	 5	 0.0002	 0.002	 0.0	 100	 100	 100	 0.0	 0.0	 1	 -360	 360;
];
