#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <map>
#include <string.h>
#include <time.h>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;


const double pi = 2.0 * acos(0);




class TangFlux{
public:
  
  double pmax = 100000; // GeV
  
  vector <double> tbl_energy;
  vector <double> tbl_range;

  TangFlux(){
    ifstream fin("standardrock.txt");
    if (fin.eof()) exit(-1);
    double ene,ran;
    while(fin >> ene >> ran){
      tbl_energy.push_back(ene);
      tbl_range.push_back(ran);
    }
  };

  double ReadTable(double dl){
    for (int i = 0; i <= tbl_range.size() - 2; i++){
      if (tbl_range[i] <= dl && dl < tbl_range[i+1]){
	double logr = log(dl);
	double logrmax = log(tbl_range[i+1]);
	double logrmin = log(tbl_range[i]);
	double logemax = log(tbl_energy[i+1]);
	double logemin = log(tbl_energy[i]);
	double loge = (logemin * (logrmax - logr) + logemax * (logr - logrmin)) / (logrmax - logrmin);
	return exp(loge);
      }
    }
    return -1;
  }

  double GetFlux(double p, double theta){

    double costheta = cos(theta);
    const double p1 = 0.102573;
    const double p2 = -0.0068287;
    const double p3 = 0.958633;
    const double p4 = 0.0407253;
    const double p5 = 0.817285;
    double costheta2 = sqrt((costheta*costheta+p1*p1+p2*pow(costheta,p3)+p4*pow(costheta,p5))/
			    (1.0+p1*p1+p2+p4));

    double p0 = (p < 1.0/costheta2 ? (3.0*p+7.0/costheta2)/10.0 : p);
    
    const double rc = 1.0e-4;
    double delta = 2.06e-3 * (950.0/costheta2 - 90.0);
    double pprod = p0 + delta;
    double a = 1.1 * pow(90.0*sqrt(costheta+0.001)/1030.0, 4.5/(pprod*costheta2));
    
    double flux;
    flux = a * 0.14 * pow(p0, -2.7) * (1.0/(1.0+1.1*pprod*costheta2/115.0)+
                                       0.054/(1.0+1.1*pprod*costheta2/850.0)+rc);

    return flux; 
  };

  double GetFlux(double p, double theta, double height){
    double h0 = 3400.0 + 1100.0 * p * cos(theta);
    double factor = exp(height/h0);
    return GetFlux(p, theta) * factor; // cm-2 sec-1 sr-1 gev-1
  };

  double GetIntegratedFlux(double pmin, double theta, double height = 0){
    double sum = 0;
    const double fact = 1.001;
    for (double p = pmin; p < pmax; p *= fact){
      double dndp = GetFlux(p, theta, height);
      sum += dndp * (p * (fact - 1.00));
    }
    return sum; // cm-2 sec-1 sr-1
  };

  double GetPenetratedFlux(double dl, double theta, double height = 0){ // dl: m.w.e
    // calculate cutoff energy
    double ecut = ReadTable(dl); // GeV
    if (ecut < 0) return 0; // exception 
    return GetIntegratedFlux(ecut, theta, height); // cm-2 sec-1 sr-1
  };

};



// This program returns the muon flux after passing through material with a given thickness.
// The opacity should be given in the unit of meter water equivallent.
// The zenith angle should be given in the unit of radian (0: vertical, pi/2: horizontal)

// The energy spectrum model used in this work is proposed by
// https://arxiv.org/abs/hep-ph/0604078 and it is used in
// several muon radiography works such as https://www.nature.com/articles/s41598-019-43527-6
// https://www.nature.com/articles/s41598-018-21423-9

// This program is for sure valid within zenith < 80 degree, opacity < 2 km.w.e.
// See Figure 3 of https://www.nature.com/articles/s41598-019-43527-6

int main(int argc, char **argv){

  // Class of Flux Model
  TangFlux rf;

  // Read Calculation Configurations
  float detsize = 1.0;
  float exposure_time = 1.0e-7; // sec 
  float numtanx = 20;
  float tanxmin = -1;
  float tanxmax = +1;
  float numtany = 10;
  float tanymin = 0;
  float tanymax = +1;
  float density_water = 1;  // g/cm3
  float height_water = 300; // meter

  cerr << " situation " << endl;
  cerr << "                   xxxxxxxxxxxx" << endl;
  cerr << "                   xxxxxxxxxxxx" << endl;
  cerr << "         ----------xx water xxx" << endl;
  cerr << "         |detector|xxxxxxxxxxxx" << endl;
  cerr << "         ----------xxxxxxxxxxxx" << endl;
  cerr << "*******************************" << endl;
   

  cerr << "Detector Size in m^2 ?: ";
  cin >> detsize;
  cerr << "Exposure Time in sec ?: ";
  cin >> exposure_time;
  cerr << "Number of Bins in X (azimuth) direction ?: ";
  cin >> numtanx;
  cerr << "Tangent(X) Min? :";
  cin >> tanxmin;
  cerr << "Tangent(X) Max? :";
  cin >> tanxmax;
  cerr << "Number of Bins in Y (elevation) direction ?: ";
  cin >> numtany;
  cerr << "Tangent(Y) Min? :";
  cin >> tanymin;
  cerr << "Tangent(Y) Max? :";
  cin >> tanymax;
  cerr << "Density of Water in g/cm3 ? :";
  cin >> density_water;
  cerr << "Height of Water in meter ? :";
  cin >> height_water;

  // Loop for each bin in angular space
  for (int ix = 0; ix < numtanx; ix++){
    float bin_tanxmin = tanxmin + (tanxmax - tanxmin) * float(ix)   / float(numtanx);
    float bin_tanxmax = tanxmin + (tanxmax - tanxmin) * float(ix+1) / float(numtanx);
    float bin_tanxcenter = (bin_tanxmin + bin_tanxmax) * 0.5;

    for (int iy = 0; iy < numtany ; iy++){
      float bin_tanymin = tanymin + (tanymax - tanymin) * float(iy)   / float(numtany);
      float bin_tanymax = tanymin + (tanymax - tanymin) * float(iy+1) / float(numtany);
      float bin_tanycenter = (bin_tanymin + bin_tanymax) * 0.5;

      // Length of muon path in the water
      float length_water = height_water * sqrt(1.0 + bin_tanxcenter * bin_tanxcenter + bin_tanycenter * bin_tanycenter); // meter
      float density_length_water = length_water * density_water; // meter water equivallent

      // Zenith Angle of muon path
      float zenith_angle = acos(bin_tanycenter / sqrt(1.0 + bin_tanxcenter * bin_tanxcenter + bin_tanycenter * bin_tanycenter));

      // calc muon flux
      float bin_muonflux = rf.GetPenetratedFlux(density_length_water, zenith_angle);

      // convert flux -> number muons
      float solid_angle = (bin_tanymax - bin_tanymin) * (bin_tanxmax - bin_tanxmin); // steradian (rough estimate. Not exact)
      float eff_area = detsize / sqrt(1.0 + bin_tanxcenter * bin_tanxcenter + bin_tanycenter * bin_tanycenter) * 10000; // m2 -> cm2
      float num_muon = bin_muonflux * solid_angle * eff_area * exposure_time;

      // Output
      cout << ix << " " 
           << iy << " " 
           << bin_tanxmin << " " 
           << bin_tanxmax << " "
           << bin_tanymin << " " 
           << bin_tanymax << " " 
           << length_water << " "
           << density_length_water << " "
           << bin_muonflux << " "
           << num_muon << endl; 
    }
  }
}






