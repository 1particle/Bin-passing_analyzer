/*
 * 	PMF_calculation.h
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

#ifndef PMF_CALCULATION_H_
#define PMF_CALCULATION_H_

class bin {
public:
  double xi, xf, W, W_sq, exp_half_W, exp_W, W_temp, den,
  sd_mean_W, sd_mean_exp_half_W, exp_2W, sd_mean_exp_W, v, v_sd, D;
  int w_count, ts_count;
  bin(double x1, double x2) : xi(x1), xf(x2), W(0), W_sq(0), exp_half_W(0), exp_W(0),
		  W_temp(0), den(0),sd_mean_W(0), sd_mean_exp_half_W(0), exp_2W(0), sd_mean_exp_W(0),
       v(0), v_sd(0), D(0),  w_count(0), ts_count(0)  { }
};

const double kB = 0.0019872041; // Boltzmann constant in kcal/(mol.K)

long file_size ();
int block_sizes(int *bs);
void calculate_time_step(vector<bin> & b, double x1, double xx1, double x2, double xx2 , double f , double ff, double bin_size, int c, int block_size, double beta);
void Work_bin_analysis (vector<bin> & b, int c, double dt);
void calculate_PMF_FR (vector<bin> & f, vector<bin> & r,  vector<bin> & w , int c, int block_size, double bin_size, double beta, double dt);
void calculate_PMF_BDFDT (vector<bin> & f, vector<bin> & r,  vector<bin> & w , double beta,  int c, int block_size, double dt);
void read_input ( vector < vector<bin> > & f,  vector < vector<bin> > & r, int c, int n, int *block_size, double bin_size, double beta);
void bin_2d_vector_creator (int c, int n, double xi, double bin_size, vector < vector<bin> > & b );
void calculate_PMF_Jarzynski (vector<bin> & f, vector<bin> & r,  vector<bin> & w , double beta,  int c, int block_size, double dt);

#endif /* PMF_CALCULATION_H_ */
