This is a parameter file for a 2-layer model with free-surface 
2 1 6 1 32.0 0.008 20 0.1
2: two domains; 1: one source; 6: six elements per wavelength; 1: ?
32.0: lengh of time T; 0.008: dt; 20: appended twenty elements on right and left edges for damping; 0.1: coefficent for damping

2.0 1.0 0.0 3.0 0.5 0.5
2.0: ?; 1.0: central frequency of ricker; 0.0: mininum frequency; 3.0: max frequency; 
0.5: length of element for minium frequency; 0.5: length of element for max frequency
The number of frequenceis is determined by: df = 1/T; In this example, T is abit larger than 32s, T=4096 * 0.008; 4096=2^12

2 -100.0 840.0 3 5.0 90.0 0.0
2: number of sources (the max can not exceed nfs); -100: x of source; 840: y of source; 3: explorational source, type of source
5.0, 90.0, 0.0: parameter for double couple source.

161 (number of receivers)
init (output file of all frequencies)
initu.dat (output file of ux)
initw.dat (output file of uy)
domain 1 2.3 5.2 3.0 27
1: NO. of domain; 2.3: pho; 5.2: vp; 3.0: vs; 27: number of coords bellow

1 23 24 27 
1 23: up boundary of domain 1; 24 27: bottom boundary of domain 1; 

1     80.0   0.00 0
2     40.0   0.00 0
3     39.5   0.0308 0
4     39.0   0.1224 0
5     38.5   0.2725 0
6     38.0   0.4775 0
7     37.5   0.7322 0
8     37.0   1.0305 0
9     36.5   1.3650 0
10    36.0   1.7275 0
11    35.5   2.1089 0
12    35.0   2.5000 0
13    34.5   2.8911 0
14    34.0   3.2725 0
15    33.5   3.6350 0
16    33.0   3.9695 0
17    32.5   4.2678 0
18    32.0   4.5225 0
19    31.5   4.7275 0
20    31.0   4.8776 0
21    30.5   4.9692 0
22    30.0   5.00 0
23   -60.0   5.00 0
24   -60.0   30.0 0
25    30.0   30.0 0
26    40.0   30.0 0
27   100.0   30.0 0
domain 2 2.5 7.6 5.0 2
1 2 3 4
1 2: up boundary of domain 2 (which should be consist with bottom boundary of domain 1 and is ignore in calculation); 3 4: bottom boundary of domain 2

1     60.0   0.0 0
2    -10.0   0.0 0
3    -10.0   30.0 0
4    -09.0   30.0 0
cccccccccccccccccccccccccccccccccccccc