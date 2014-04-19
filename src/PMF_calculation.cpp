/*
 * 	PMF_calculation.cpp
 * 	Copyright 2014 Mostafa Nategholeslam

	This file is part of Bin-passing_analyzer.

    Bin-passing_analyzer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Bin-passing_analyzer is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
using namespace std;
#include <time.h>
#include <iomanip>
#include <locale>
#include <sstream>
#include <sys/stat.h>
//#include <sys/types.h>
#include <cstring>
#include "PMF_calculation.h"



long file_size ()  {
  long begin,end;
  ifstream input ("FR.txt");
  begin = input.tellg();
  input.seekg (0, ios::end);
  end = input.tellg();
  input.close();
  return end-begin;
}

void calculate_time_step(vector<bin> & b, double x1, double xx1, double x2, double xx2 , double f , double ff, double bin_size, int c, int block_size, double beta)
{
  // Positive value for variable 'bin_size' means forward direction is positive (from smaller to larger values of x) and negative means the reverse.
  double X1 = x1 - xx1;
      double X2 = x2 - xx2;
      double dx= (x2-x1);
      double dxx= (xx2-xx1);
      double dX = X2 - X1;
      int i= (int) floor ((X1 - b[0].xi)/bin_size);

      if (i>=0 && i<(int)b.size())
      if ( (bin_size > 0 && b[i].xi < X1 && X1  < b[i].xf &&   b[i].xi < X2  &&  X2  < b[i].xf )
            || (bin_size < 0 && b[i].xi > X1 && X1  > b[i].xf &&   b[i].xi > X2  &&  X2  > b[i].xf ))
            {
                    b[i].ts_count ++;
                    b[i].W_temp += (f*dx + ff*dxx);
                    b[i].den += fabs((dX)/bin_size);
                    b[i].v += fabs(dX);
                    b[i].v_sd += dX*dX;
                    if(b[i].den >= block_size)         {
                            double W = b[i].W_temp/b[i].den;
                            b[i].W += W;
                            b[i].W_sq += W*W;
                            b[i].exp_half_W += exp(-0.5000*beta*W);
                            b[i].exp_W += exp(- beta*W) ;
                            b[i].exp_2W += exp(- 2*beta*W) ;
                            b[i].w_count ++;

                        b[i].den =0;
                        b[i].W_temp=0;
                    }
                          }
}


void Work_bin_analysis (vector<bin> & b, int c, double dt)
{
  for (int i=0; i < c ; i++)             {
      if ( b[i].w_count &&  b[i].w_count )    {
            b[i].v /= (double)((double)b[i].ts_count * dt);
            b[i].v_sd /= (double)((double)b[i].ts_count * dt*dt);
            b[i].v_sd = sqrt(b[i].v_sd-b[i].v*b[i].v)/b[i].ts_count;
            b[i].W /= (double)b[i].w_count;
            b[i].W_sq /= (double)b[i].w_count;
            b[i].exp_half_W /= (double)b[i].w_count;
            b[i].exp_W /= (double)b[i].w_count;
            b[i].exp_2W /= (double)b[i].w_count;
            b[i].sd_mean_W = sqrt((b[i].W_sq -b[i].W*b[i].W)/(double)b[i].w_count);
            b[i].sd_mean_exp_half_W = sqrt((b[i].exp_W -b[i].exp_half_W*b[i].exp_half_W)/(double)b[i].w_count);
            b[i].sd_mean_exp_W = sqrt((b[i].exp_2W -b[i].exp_W*b[i].exp_W)/(double)b[i].w_count);

  }
  }
}


void calculate_PMF_FR (vector<bin> & f, vector<bin> & r,  vector<bin> & w , int c, int block_size, double bin_size, double beta, double dt)
{
  ofstream results;
  results.precision(18);
  ostringstream convert;
  convert << block_size;
 string filename = "PMFs_NEW/PMF_FR_block_size_"+ convert.str() + ".txt";
  results.open(filename.c_str());
  for (int i=0; i < c ; i++)             {

            w[i].W = (f[i].W-r[i].W)/2;
            w[i].D = (f[i].v+r[i].v) * bin_size / ((f[i].W+r[i].W) * beta);
            w[i].D /= 100.0000;   // To convert from (A^2)/ns to units of (10^-5 cm^2)/s
            w[i].sd_mean_W = sqrt (pow(f[i].sd_mean_W,2)+pow(r[i].sd_mean_W,2));

            if (i>0)        {
              w[i].W += w[i-1].W;
              w[i].sd_mean_W = sqrt(pow(w[i].sd_mean_W,2)+pow(w[i-1].sd_mean_W,2));
            }
  }
            // Getting the errors from the opposite size and combining the errors with appropriate weights
            for (int i=c-1; i >= 0 ; i --)         {
              if ( f[i].w_count && r[i].w_count )        {
                if (i<c-1)	{
                	w[i].W_temp = sqrt (pow(f[i+1].sd_mean_W,2)+pow(r[i+1].sd_mean_W,2));
                	w[i].W_temp = sqrt(pow(w[i].W_temp,2) + pow(w[i+1].W_temp,2));
                }
                else
                	w[i].W_temp = 0;

                 w[i].sd_mean_W = (w[i].sd_mean_W * w[i].W_temp)/sqrt(pow(w[i].sd_mean_W, 2)+ pow(w[i].W_temp,2));

              }
            }

            results << "Block_size= " << block_size << endl <<endl;
            results << "x_f_bin   PMF_FR  d(PMF_FR)  D(z) v_F  d(v_F)  v_R  d(v_R)  W_d   "
            		"w_F   d(w_F)  w_R  d(w_R)   "
            		"w_R_block_counts   w_R_block_counts  t_forward(ns  t_reverse(ns)" <<endl;
            results<< w[0].xi << "  "<<"0.000" <<"   "<<"0.000"<<endl;
            for (int i=0; i < c ; i++)  {
				 results<< w[i].xf << "  "<<w[i].W<<"   "<<w[i].sd_mean_W<< "   "<< w[i].D <<"   "
						 << f[i].v<<"  "<<f[i].v_sd<<"  "<< r[i].v <<"    "<<r[i].v_sd<<"   "<< (f[i].W+r[i].W)/2<<"   "
						 <<f[i].W<<"  "<<f[i].sd_mean_W<<"  "<<r[i].W<<"   "<<r[i].sd_mean_W<<"  "
						 <<"   "<<f[i].w_count<<"  "<<r[i].w_count <<"  "
						 <<(double)f[i].ts_count*dt<<"  "<<(double)r[i].ts_count*dt<<endl;
                                 }
            results.close();

}


void calculate_PMF_BDFDT (vector<bin> & f, vector<bin> & r,  vector<bin> & w , double beta,  int c, int block_size, double dt)
{

  ofstream results;
  results.precision(18);
    ostringstream convert;
    convert << block_size;
   string filename = "PMFs_NEW/PMF_BD-FDT_block_size_"+ convert.str() + ".txt";
    results.open(filename.c_str());

  for (int i=0; i < c ; i++)             {

      w[i].W= - log(f[i].exp_half_W/r[i].exp_half_W)/beta;
      w[i].sd_mean_W = sqrt(pow(f[i].sd_mean_exp_half_W/f[i].exp_half_W,2)  +   pow(r[i].sd_mean_exp_half_W/r[i].exp_half_W,2))/beta;
          if (i>0)        {
              w[i].W += w[i-1].W;
              w[i].sd_mean_W = sqrt(pow(w[i].sd_mean_W,2)+pow(w[i-1].sd_mean_W,2));
          }    ////////////////////////////////////////////////////////////////////////////////////////////

  }
            // Getting the errors from the opposite size and combining the errors with appropriate weights
            for (int i=c-1; i >= 0 ; i --)         {
              if ( f[i].w_count && r[i].w_count )        {
                  if (i<c-1)	{
                  	w[i].W_temp = sqrt (pow(f[i+1].sd_mean_W,2)+pow(r[i+1].sd_mean_W,2));
                  	w[i].W_temp = sqrt(pow(w[i].W_temp,2) + pow(w[i+1].W_temp,2));
                  }
                  else
                  	w[i].W_temp = 0;
				   w[i].sd_mean_W = (w[i].sd_mean_W * w[i].W_temp)/sqrt(pow(w[i].sd_mean_W, 2)+ pow(w[i].W_temp,2));

              }
            }

            results << "Block_size= " << block_size << endl <<endl;
            results << "x_f_bin   PMF_BD-FDT d(PMF_BD-FDT)  <e^(-beta*w_F)>*<e^(-beta*w_R)>   w_F_block_counts  w_R_block_counts  t_forward(ns) t_reverse(ns)" <<endl;
            results<< w[0].xi << "  "<<"0.000" <<"   "<<"0.000"<<endl;
                     for (int i=0; i < c ; i++) {
                     results<< w[i].xf << "  "<<w[i].W<<"   "<<w[i].sd_mean_W<< "  "<<f[i].exp_W*r[i].exp_W<< "   ";
                     results<< f[i].w_count<<"  "<<r[i].w_count<<"  "<<(double)f[i].ts_count*dt<<"  "<<(double)r[i].ts_count*dt<<endl;
                     }
            results.close();

}

void calculate_PMF_Jarzynski (vector<bin> & f, vector<bin> & r,  vector<bin> & w , double beta,  int c, int block_size, double dt)
{

  ofstream results;
  results.precision(18);
    ostringstream convert;
    convert << block_size;
   string filename = "PMFs_NEW/PMF_Jarzynski_block_size_"+ convert.str() + ".txt";
    results.open(filename.c_str());

  for (int i=0; i < c ; i++)             {

      w[i].W= - log(f[i].exp_W)/beta;
      w[i].sd_mean_W = (f[i].sd_mean_exp_W/f[i].exp_W)/beta;
          if (i>0)        {
              w[i].W += w[i-1].W;
              w[i].sd_mean_W = sqrt(pow(w[i].sd_mean_W,2)+pow(w[i-1].sd_mean_W,2));
          }    ////////////////////////////////////////////////////////////////////////////////////////////

  }
            // Getting the errors from the opposite size and combining the errors with appropriate weights
            for (int i=c-1; i >= 0 ; i --)         {
              if ( f[i].w_count && r[i].w_count )        {
                  if (i<c-1)	{
                  	w[i].W_temp = sqrt (pow(f[i+1].sd_mean_W,2)+pow(r[i+1].sd_mean_W,2));
                  	w[i].W_temp = sqrt(pow(w[i].W_temp,2) + pow(w[i+1].W_temp,2));
                  }
                  else
                  	w[i].W_temp = 0;
                  w[i].sd_mean_W = (w[i].sd_mean_W * w[i].W_temp)/sqrt(pow(w[i].sd_mean_W, 2)+ pow(w[i].W_temp,2));

              }
            }

            results << "Block_size= " << block_size << endl <<endl;
            results << "x_f_bin   PMF_Jarzynski d(PMF_Jarzynski)  <e^(-beta*w_F)>*<e^(-beta*w_R)>   w_F_block_counts  w_R_block_counts  t_forward(ns) t_reverse(ns)"<<endl;
            results<< w[0].xi << "  "<<"0.000" <<"   "<<"0.000"<<endl;
            for (int i=0; i < c ; i++)  {
                     results<< w[i].xf << "  "<<w[i].W<<"   "<<w[i].sd_mean_W<< "  "<<f[i].exp_W*r[i].exp_W<< "   ";
                     results<< f[i].w_count<<"  "<<r[i].w_count<<"  "<<(double)f[i].ts_count*dt<<"  "<<(double)r[i].ts_count*dt<<endl;
            }
                     results.close();

}


void read_input ( vector < vector<bin> > & f,  vector < vector<bin> > & r, int c, int n, int *block_size, double bin_size, double beta)              {

  long fileSize = file_size ();
  long current=1;
  int percent = 0;
  cout <<"Reading values from FR.txt .... \n\n\n";
  string line;
  char* st;
  ifstream FR;
  FR.open("FR.txt");
  FR.precision(15);
  while ( ! FR.eof() )
          {
                  getline (FR,line);
                  current += (int)line.size()+1;  //Need to add 1 for the newline character which is not included in the line character. (?????)
                  if ( ((current*100)/fileSize) > percent)       {
                      percent ++;
                      cout <<"\r"<<percent<<"% read..... "<< std::flush;
                      }
                  double x1 , x2, F,  xx1, xx2, FF;
                  x1= strtod(line.c_str(),&st);
                  x2= strtod(st,&st);
                  F= strtod(st,&st);
                  xx1= strtod(st,&st);
                  xx2= strtod(st,&st);
                  FF= strtod(st,&st);

                  double dx= (x2 - xx2 - x1 + xx1);
                  if (dx*bin_size > 0)
                    for (int a=0; a<n;a++)
                    calculate_time_step(f[a],  x1,  xx1,  x2,  xx2 , F, FF ,  bin_size, c, block_size[a] , beta);
                  else
                    for (int a=0; a<n;a++)
                    calculate_time_step(r[a],  x1,  xx1,  x2,  xx2 , F, FF ,  bin_size, c, block_size[a], beta);

          }

  FR.close();

}

int block_sizes(int *bs)      {
  int n=1;
  cout <<"Enter the desired number of different block sizes (maximum of 20):  ";
  cin >> n;
   for (int a=0; a<n; a++)      {
      cout <<a+1<<"  Enter the desired block size:  ";
      cin >> bs[a];
  }
  return n;
}

void bin_2d_vector_creator (int c, int n, double xi, double bin_size, vector < vector<bin> > & b )    {
for (int i=0; i<n; i++)   // n is the number of different block sizes, and c is the number of bins along the reaction path
      {
          vector<bin> row;
          for (int j=0; j<c; j++) {
                          bin B(xi + j*bin_size, xi + (j+1)*bin_size);
                          row.push_back(B);
                  }
          b.push_back(row);
}
}
