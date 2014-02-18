/*
 * 	Bin-passing_analyzer.cpp
 * 	Copyright 2014  © Mostafa Nategholeslam

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

int main(int argc, char* argv[]) {

//  gtk_init(&argc, &argv);
// Win();
//     gtk_main();

        mkdir("PMFs_NEW",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        double init, fin, T=310, dt=2.0000;
        int c;
        cout<<endl<<"Enter the initial value of the range of reaction coordinate:  ";
        cin>>init;
        cout<<"Enter the final value of the range of reaction coordinate:  ";
        cin>>fin;
        cout<<"How many bins do you want along this range of reaction coordinate?  ";
        cin>>c;
        cout <<"At what temperature (in Kelvin) did you run this simulation? ";
        cin >> T;
        cout <<"Enter the length of the time-step used for this simulation (in femtoseconds): ";
        cin >> dt;
        dt *= 1.00000e-6;

        double bin_size=(fin-init)/(double)c;   //ATTENTION: This can be negative.
        cout <<"\n Bin size = "<<bin_size<<endl<<endl;
        int block_size[20]={0};
        int n = block_sizes(block_size);

        vector <vector<bin> > f, r, W_FR, W_BDFDT,W_Jarzynski;
        bin_2d_vector_creator (c, n, init, bin_size, f );
        bin_2d_vector_creator (c, n, init, bin_size, r);
        bin_2d_vector_creator (c, n, init, bin_size, W_FR );
        bin_2d_vector_creator (c, n, init, bin_size, W_BDFDT );
        bin_2d_vector_creator (c, n, init, bin_size, W_Jarzynski );


        time_t seconds;
         seconds = time (NULL);

         double beta = 1/(kB*T);

        cout.precision(18);

        read_input (f,  r,  c, n, block_size, bin_size, beta);
        for (int i=0; i<n; i++)         {
        Work_bin_analysis(f[i],c, dt);
        Work_bin_analysis(r[i],c, dt);
        calculate_PMF_FR ( f[i],  r[i],  W_FR[i] , c, block_size[i], fabs(bin_size), beta, dt);
        calculate_PMF_BDFDT (f[i],  r[i],  W_BDFDT[i] , beta, c, block_size[i], dt);
        calculate_PMF_Jarzynski (f[i],  r[i],  W_Jarzynski[i] , beta, c, block_size[i], dt);
        }


        seconds = time (NULL) - seconds;
        cout<<endl<<endl<<"Done!"<<endl<<"And it took "<<seconds<<" seconds to do it!"<<endl<<endl;


        return 0;
}

