# calc-muflux-dam
a simple code to simualte the muon flux using r-nishi-eri/calc-muflux-tang

Usage
1) Compile using g++
g++ calc-muflux-dam.cpp -o calc-muflux-dam

2) Execute the executable
./calc-muflux-dam

3) Input the values
 situation
                   xxxxxxxxxxxx
                   xxxxxxxxxxxx
         ----------xx water xxx
         |detector|xxxxxxxxxxxx
         ----------xxxxxxxxxxxx
*******************************
Detector Size in m^2 ?: 1
Exposure Time in sec ?: 1000000
Number of Bins in X (azimuth) direction ?: 20
Tangent(X) Min? :-1
Tangent(X) Max? :1
Number of Bins in Y (elevation) direction ?: 10
Tangent(Y) Min? :0
Tangent(Y) Max? :1
Density of Water in g/cm3 ? :1
Height of Water in meter ? :200

4) You get the output
0 0 -1 -0.9 0 0.1 276.043 276.043 2.02144e-05 1464.58
0 1 -1 -0.9 0.1 0.2 277.489 277.489 2.86474e-05 2064.76
0 2 -1 -0.9 0.2 0.3 280.357 280.357 3.35542e-05 2393.68
0 3 -1 -0.9 0.3 0.4 284.605 284.605 3.55966e-05 2501.47
0 4 -1 -0.9 0.4 0.5 290.172 290.172 3.57372e-05 2463.17
0 5 -1 -0.9 0.5 0.6 296.985 296.985 3.47103e-05 2337.51
0 6 -1 -0.9 0.6 0.7 304.959 304.959 3.29924e-05 2163.72
0 7 -1 -0.9 0.7 0.8 314.006 314.006 3.09118e-05 1968.87 ...

$1 index x
$2 index y 
$3 Min. of tan_x
$4 Max. of tan_x
$5 Min. of tan_y
$6 Max. of tan_y
$7 Length of water (m)
$8 Density Length (meter water equivallent)
$9 Muon Flux (cm^-2 sec^-1 steradian^-1)
$10 Number of Muons
