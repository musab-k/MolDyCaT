#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;


double phonon_dos(string seed, const vector<double>& freqs,const double min_w, const double max_w, const double w_step, const double smoothing, int branches)
{
    ofstream of;
    of.open(string(seed + ".dos").c_str());

    double w_range = max_w - min_w;
    int n = floor(w_range/w_step);
    for(int i = 0; i <= n; i++)
    {
        double w = min_w + i * w_step;
        double g = 0.0;

        double pf = 1.0/(smoothing * sqrt(2.0*M_PI));

        for(unsigned int s = 0; s < freqs.size(); s++)
        {
                double w0 = freqs[s];
                g += pf * exp(-1.0*(w-w0)*(w-w0)/(smoothing*smoothing));
        }
        g /= (8.0*M_PI*M_PI*M_PI * branches);

        of << w << '\t' << g << endl;
    }

    of.close();
    return 0;
}




int main(int argc, char* argv[])
{
    if(argc == 1)
    {
        cerr << "No seed specified" << endl; 
        return 1;
    }
   string seed(argv[1]);

  ifstream in;
  string paramfile = seed + ".freq";
  in.open(paramfile.c_str());
  if(!in){
    cerr << "Input file could not be found \n";
    return 1;
  }
  string str;
  vector<double> freqs(0);

  string q = "q-point";
  int branches = 0;

  while(!in.eof())
  {
    getline(in,str);
    if (str.compare(0, q.length(), q) == 0)
    {
        branches++;
        continue;
    }
    else
    {
        if(str.compare(0,1," ") != 0 && str.compare(0,1,"\n") != 0 )
            freqs.push_back(atof(str.c_str()));
    }
      
  }
  double min_w, max_w, w_step, smoothing;
  cout << "Please enter lowest frequency" << endl;
  cin >> min_w;
  cout << "Please enter highest frequency" << endl;
  cin >> max_w;
  cout << "Please enter frequency step" << endl;
  cin >> w_step;
  cout << "Please enter smoothing factor" << endl;
  cin >> smoothing;
 in.close();

 phonon_dos(seed,freqs,min_w,max_w,w_step,smoothing,branches);
 cout << "D.O.S written to " << seed << ".dos" << endl;
}
/*

double hbar = 1.0;

double phonon_free_energy(const vec freqs,const unitcell& cell, const double T)
{
    const double beta = 1.0 / (kB * T);

    double boltz_sum = 0.0;
    for(int i = 0; i < freqs.count(); i++)
    {
        boltz_sum += exp(-beta*hbar*freqs[i]);
    }
    return -ln(z)/(beta*cell.volume());
}

double phonon_mode_cv(const double w, const double beta, unitcell& cell)
{
    double bf = exp(beta*hbar*w);
    return  (beta/T) *hbar*hbar*w*w / ((bf-1.0)*(bf-1.0) * cell.volume());
}


double phonon_heat_capacity(const vec freqs,const unitcell& cell, const double T)
{
    const double beta = 1.0 / (kB * T);

    double cv = 0.0;
    for(int i = 0; i < freqs.count(); i++)
    {
       cv += phonon_mode_cv(freqs[i],cell,beta);
    }

    return cv;

}
*/
