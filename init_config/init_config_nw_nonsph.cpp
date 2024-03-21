//
//  Created by Tsai-Wei Lin on 08/04/21.
//

#include <iostream>    //header for basic io
#include <cmath>       //header for math functions
#include <fstream>     //header for file io
#include <cstdlib>     //header for srand
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;

pid_t getpid(void);
#define PI 3.14159265358979323846
#define IM1 2147483563 /* what does this number mean? */
#define IM2 2147483399 /*this is part of the ran numb gen*/
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
double ran1(long *idum);
double gasdev(long *idum);
long initRan();

int main(int argc, const char * argv[]) {

    int N_chain, N_mono, N_cross, N_dangle, Nrod, AR;
    int ranB, id_chosen, id_prev, id_new, id_chosen_all, id_new_all, id_prev_all;
    int testflag = 0, flag2 = 0;
    double f_cross_b, f_cross, box_length, half_box, avg_bondlen, len_pene, Dp;
    double theta, phi, vecx, vecy, vecz, dx, dy, dz;
    char null[1024], fout_conf[1024], fn_seed[1024];
    long *idum = (long *) malloc(sizeof(long));
    long seed_num;
    ofstream out_config, out_seedls; 

    /// Command line inputs

    N_chain = atoi(argv[1]);    // number of chains
    N_mono = atoi(argv[2]);     // number of beads in each chain  
    N_cross = atoi(argv[3]);    // number of reactive bead on each chain
    N_dangle = atoi(argv[4]);   // number of dangling bead on each reative bead
    Nrod = atoi(argv[5]);       // number of rod-like penetrants
    AR = atoi(argv[6]);         // how many bead a penetrants is composed of
    len_pene = atof(argv[7]);   // chain length between penetrant beads
    Dp = atof(argv[8]);         // penetrant bead diameter

    /// Vectors - Polymer chains
    vector<vector<double> > chainrx(N_chain);
    vector<vector<double> > chainry(N_chain);
    vector<vector<double> > chainrz(N_chain);
    vector<vector<double> > charge(N_chain);
    vector<vector<int> > atype(N_chain);
    vector<vector<int> > btype(N_chain);
    vector<vector<int> > bfirst(N_chain);
    vector<vector<int> > bsecond(N_chain);
    vector<vector<int> > angtype(N_chain);
    vector<vector<int> > angfirst(N_chain);
    vector<vector<int> > angsecond(N_chain);
    vector<vector<int> > angthird(N_chain);

    /// Vectors - Rod penetrant
    vector<vector<double> > penerx(Nrod);
    vector<vector<double> > penery(Nrod);
    vector<vector<double> > penerz(Nrod);
    vector<vector<double> > pene_q(Nrod);
    vector<vector<int> > patype(Nrod);
    vector<vector<int> > pbtype(Nrod);
    vector<vector<int> > pbfirst(Nrod);
    vector<vector<int> > pbsecond(Nrod);
    vector<vector<int> > pangtype(Nrod);
    vector<vector<int> > pangfirst(Nrod);
    vector<vector<int> > pangsecond(Nrod);
    vector<vector<int> > pangthird(Nrod);


    f_cross_b = (double) N_cross/(N_mono-N_cross);     // definition before crosslinking
    f_cross = (double) N_cross/(N_mono-2*N_cross);     // definition after crosslinking
    box_length = pow(((N_mono+N_cross*N_dangle)*N_chain/0.85), 1.0/3.0);
    avg_bondlen = 1.3000;
    half_box = box_length/2.0;

    time_t theTime;
    time(&theTime);
    /// Seeds the pseudo-random number generator used by rand() with the value seed.
    srand(time(0));
    /// Initialize Random Variable
    *idum = (-1)*initRan();
    seed_num = (-1)*(*idum);
    //printf("SEED = %ld\n", seed_num);

    sprintf(fout_conf, "./f_%.3lf/initcfg-Np%d-AR%d-Dp%.3lf-lp%.3f/f_%.3lf-%ld-Nm%d-Nc%d.txt", f_cross, Nrod, AR, Dp, len_pene, f_cross, seed_num, N_mono, N_chain);
    sprintf(fn_seed, "./f_%.3lf/initcfg-Np%d-AR%d-Dp%.3lf-lp%.3f/Seed_list.txt", f_cross, Nrod, AR, Dp, len_pene);
    
    /// Save the seed
    out_seedls.open(fn_seed, ios::out | ios::app);
    out_seedls << seed_num << endl;


    
    /// Bead positions - polymer chains
    for (int i=0; i<N_chain; ++i){
        for (int j=0; j<N_mono; ++j){
            testflag = 0;
            
            while(testflag==0){
                testflag = 1;

                if (j==0){
                    vecx = box_length*(ran1(idum))-half_box;
                    vecy = box_length*(ran1(idum))-half_box;
                    vecz = box_length*(ran1(idum))-half_box;
                    // vecx = vecy = vecz = 1.2;
                    chainrx[i].push_back(vecx);
                    chainry[i].push_back(vecy);
                    chainrz[i].push_back(vecz);
                    atype[i].push_back(1);
                    charge[i].push_back(0.0);
                }
                else{
                    testflag = 1;
                    flag2 = 1;
                    theta = ran1(idum)*2.0*3.14158; //theta ranges from 0 to 2pi
                    phi = acos(2.0*ran1(idum)-1.0); //phi ranges from 0 to pi
                    vecx = chainrx[i][j-1] + avg_bondlen*cos(theta)*sin(phi);
                    vecy = chainry[i][j-1] + avg_bondlen*sin(theta)*sin(phi);
                    vecz = chainrz[i][j-1] + avg_bondlen*cos(phi);
                    vecx -= box_length*round(vecx/box_length);
                    vecy -= box_length*round(vecy/box_length);
                    vecz -= box_length*round(vecz/box_length);


                    /// Prevent bead in the same chain from overlapping (Do we really need this step if we are using soft potential first?)
                    for (int k=0; k<j; ++k){
                        dx = chainrx[i][k] - vecx;
                        dy = chainry[i][k] - vecy;
                        dz = chainrz[i][k] - vecz;
                        dx -= box_length*round(dx/box_length);
                        dy -= box_length*round(dy/box_length);
                        dz -= box_length*round(dz/box_length);
                        if(dx*dx+dy*dy+dz*dz<1.0){
                            testflag = 0;
                            flag2 = 0;
                        }
                    }

                    if (flag2 ==1){
                        chainrx[i].push_back(vecx);
                        chainry[i].push_back(vecy);
                        chainrz[i].push_back(vecz);
                        atype[i].push_back(1);
                        charge[i].push_back(0.0);
                    }
                }
            }
            //cout << "Found position for atom with Nchain = " << i <<", Nmono = " << j << endl;
        }
        //cout << "done assigning back bone monomer for chain number " << i+1 << endl;
        /// Store bond info for backbone monomers
        for (int j=0; j<N_mono-1; ++j){
            btype[i].push_back(1);              //set bond between backbone monomers as bond type 1
            bfirst[i].push_back(i*(N_mono+N_cross*N_dangle)+j);         //store first atom number in this bond
            bsecond[i].push_back(i*(N_mono+N_cross*N_dangle)+j+1);      //store second atom number in this bond
        }   
        //cout << "done saving bond info for back bone monomer for chain number " << i+1 << endl;
        for (int j=0; j<N_mono-2; ++j){
            angtype[i].push_back(1);             //set angle among backbone monomers as angle type 1
            angfirst[i].push_back(i*(N_mono+N_cross*N_dangle)+j);       //store first atom number in this angle
            angsecond[i].push_back(i*(N_mono+N_cross*N_dangle)+j+1);    //store second (center) atom number in this angle
            angthird[i].push_back(i*(N_mono+N_cross*N_dangle)+j+2);     //store third atom number in this angle        
        }
        //cout << "done saving angle info for back bone monomer for chain number " << i+1 << endl;
        

        /// For each loop, we have to 
        /// 1) randomly choose backbone monomer 
        /// 2) store bond, angle, and position info for N_dangle new atoms where the last one is crosslink bead
        for(int j=0; j<N_cross; ++j){
            /// 1) Randomly choose backbone monomer
            /// 1-a) Set the chosen monomer as atom type 2
            /// 1-b) Set the new bead as atom type 3 (if the new bead is crosslink bead, its type is 4)
            testflag = 0;
            while(testflag==0){
                ranB = rand()%(N_mono-2);
                id_chosen = ranB +1;
                //cout << "the chosen monomer is number " << id_chosen << endl;
                //cout << atype[i][id_chosen] << endl;
                /// avoid choosing the same backbone monomer
                if (atype[i][id_chosen]!=2){
                    id_chosen = ranB +1;
                    testflag = 1;
                }
            }
            atype[i][id_chosen] = 2;

            /// 2) store bond, angle, and position info for N_dangle new atoms where the last one is crosslink bead

            /// Create N_dangle new beads starting from the chosen bead
            for (int k=0; k<N_dangle; ++k){
                id_prev = id_new;
                id_new = N_mono + N_dangle*j + k; // update atom numer for new-formed atom
                
                //For the bead connected to the chosen bead: +1 bond, +2 angles
                if (k==0){
                    id_prev = id_chosen;
                    atype[i].push_back(3);       //set dangling crosslink bead as type 3
                    id_chosen_all = id_chosen + i*(N_mono+N_cross*N_dangle);    // ID of the chosen backbone monomer in initial file
                    id_new_all = id_new + i*(N_mono+N_cross*N_dangle);          // ID of this new crosslink bead in initial file
                    
                    /// 2-a) Store bond info for the chosen backbone monomer and crosslink bead
                    btype[i].push_back(2);              //set bond between backbone monomer and crosslink bead as bond type 2
                    bfirst[i].push_back(id_chosen_all); //store the ID of the chosen backbone atom as first atom in the new bond
                    bsecond[i].push_back(id_new_all);   //store the ID of the new-formed crosslink bead as second atom in the new bond

                    angtype[i].push_back(2);                    //set angle between backbone monomer and crosslink bead as angle type 2
                    angfirst[i].push_back(id_chosen_all -1);    //store the ID of the backbone atom before the chosen backbone atom as first atom in the new angle
                    angsecond[i].push_back(id_chosen_all);      //store the ID of the chosen backbone atom as second (center) atom in the new angle
                    angthird[i].push_back(id_new_all);          //store the ID of the new-formed crosslink bead as third atom in the new angle

                    angtype[i].push_back(2);                    //set angle between backbone monomer and crosslink bead as angle type 2
                    angfirst[i].push_back(id_chosen_all +1);    //store the ID of the backbone atom after the chosen backbone atom as first atom in the new angle
                    angsecond[i].push_back(id_chosen_all);      //store the ID of the chosen backbone atom as second (center) atom in the new angle
                    angthird[i].push_back(id_new_all);          //store the ID of the new-formed crosslink bead as third atom in the new angle
                }

                /// For the second dangling bead: +1 bond, +1 angles
                else if (k==1){
                    
                    atype[i].push_back(3);       //set dangling crosslink bead as type 3
                    id_prev_all = id_prev + i*(N_mono+N_cross*N_dangle);
                    id_chosen_all = id_chosen + i*(N_mono+N_cross*N_dangle);    // ID of the chosen backbone monomer in initial file
                    id_new_all = id_new + i*(N_mono+N_cross*N_dangle);          // Number of this new crosslink bead in initial file
                    
                    /// 2-a) Store bond info for the chosen backbone monomer and crosslink bead
                    btype[i].push_back(2);              //set bond between first and second dangling bead as bond type 2
                    bfirst[i].push_back(id_prev_all);   //store the ID of the first dangling bead as first atom in the new bond
                    bsecond[i].push_back(id_new_all);   //store the ID of the second dangling bead as second atom in the new bond

                    angtype[i].push_back(2);                  
                    angfirst[i].push_back(id_chosen_all);    //store the ID of the chosen backbone bead as first atom in the new angle
                    angsecond[i].push_back(id_prev_all);     //store the ID of the first dangling as second (center) atom in the new angle
                    angthird[i].push_back(id_new_all);       //store the ID of the second dangling bead as third atom in the new angle
                }

                //For the last dangling bead which is crosslink bead: +1 bond, +1 angles, and set its atom type as 4
                else if (k==(N_dangle-1)){
                    atype[i].push_back(4);       //set the last dangling crosslink bead as type 4

                    id_prev_all = id_prev + i*(N_mono+N_cross*N_dangle);
                    id_new_all = id_new + i*(N_mono+N_cross*N_dangle);          // Number of this new crosslink bead in initial file
                    
                    /// 2-a) Store bond info for the chosen backbone monomer and crosslink bead
                    btype[i].push_back(2);              //set bond between first and second dangling bead as bond type 2
                    bfirst[i].push_back(id_prev_all);   //store the ID of the first dangling bead as first atom in the new bond
                    bsecond[i].push_back(id_new_all);   //store the ID of the second dangling bead as second atom in the new bond

                    angtype[i].push_back(2);                  
                    angfirst[i].push_back(id_prev_all -1);   //store the ID of the chosen backbone bead as first atom in the new angle
                    angsecond[i].push_back(id_prev_all);     //store the ID of the first dangling as second (center) atom in the new angle
                    angthird[i].push_back(id_new_all);       //store the ID of the second dangling bead as third atom in the new angle
                
                }

                else{
                    atype[i].push_back(3);       //set the dangling crosslink bead as type 3

                    id_prev_all = id_prev + i*(N_mono+N_cross*N_dangle);
                    id_new_all = id_new + i*(N_mono+N_cross*N_dangle);          // Number of this new crosslink bead in initial file
                    
                    /// 2-a) Store bond info for the chosen backbone monomer and crosslink bead
                    btype[i].push_back(2);              //set bond between first and second dangling bead as bond type 2
                    bfirst[i].push_back(id_prev_all);   //store the ID of the first dangling bead as first atom in the new bond
                    bsecond[i].push_back(id_new_all);   //store the ID of the second dangling bead as second atom in the new bond

                    angtype[i].push_back(2);                  
                    angfirst[i].push_back(id_prev_all -1);   //store the ID of the chosen backbone bead as first atom in the new angle
                    angsecond[i].push_back(id_prev_all);     //store the ID of the first dangling as second (center) atom in the new angle
                    angthird[i].push_back(id_new_all);       //store the ID of the second dangling bead as third atom in the new angle
                
                }

                /// 2-b) Create one dangling crosslink bead and store its position info
                testflag = 0;
                while(testflag==0){
                    testflag = 1;
                    theta = ran1(idum)*2.0*3.14158; //theta ranges from 0 to 2pi
                    phi = acos(2.0*ran1(idum)-1.0); //phi ranges from 0 to pi
                    vecx = chainrx[i][id_prev] + avg_bondlen*cos(theta)*sin(phi);
                    vecy = chainry[i][id_prev] + avg_bondlen*sin(theta)*sin(phi);
                    vecz = chainrz[i][id_prev] + avg_bondlen*cos(phi);
                    vecx -= box_length*round(vecx/box_length);
                    vecy -= box_length*round(vecy/box_length);
                    vecz -= box_length*round(vecz/box_length);
                    chainrx[i].push_back(vecx);
                    chainry[i].push_back(vecy);
                    chainrz[i].push_back(vecz);

                }
            }
        }
        //cout << "done chain number " << i+1 << endl;
    }
    
    // for (int i=0; i<N_chain; ++i){
    //     cout << "size of chainrx of chain "<< i+1 << "= " << chainrx[i].size() << endl;
    // }

    /// Bead positions - rod penetrants
    for (int i=0; i<Nrod; ++i){

        /// Penetrant - atoms
        for (int j=0; j<AR; ++j){
            if(j==0){
                vecx = box_length*(ran1(idum))-half_box;
                vecy = box_length*(ran1(idum))-half_box;
                vecz = box_length*(ran1(idum))-half_box;
                penerx[i].push_back(vecx);
                penery[i].push_back(vecy);
                penerz[i].push_back(vecz);
                pene_q[i].push_back(0.0);
                patype[i].push_back(6);
            }
            else{
                vecx = penerx[i][j-1] + len_pene;
                vecy = penery[i][j-1];
                vecz = penerz[i][j-1];
                vecx -= box_length*round(vecx/box_length);
                vecy -= box_length*round(vecy/box_length);
                vecz -= box_length*round(vecz/box_length);
                penerx[i].push_back(vecx);
                penery[i].push_back(vecy);
                penerz[i].push_back(vecz);
                pene_q[i].push_back(0.0);
                patype[i].push_back(6);
            }
        }

        /// Penetrant - Bonds
        for (int j=0; j<AR-1; ++j){
            pbtype[i].push_back(4);
            pbfirst[i].push_back((N_chain*(N_mono+N_cross*N_dangle)) + i*AR+j);
            pbsecond[i].push_back((N_chain*(N_mono+N_cross*N_dangle)) + i*AR+j +1);
        }

        /// Penetrant - Angles
        for (int j=0; j<AR-2; ++j){
            pangtype[i].push_back(4);
            pangfirst[i].push_back((N_chain*(N_mono+N_cross*N_dangle)) + i*AR+j);
            pangsecond[i].push_back((N_chain*(N_mono+N_cross*N_dangle)) + i*AR+j +1);
            pangthird[i].push_back((N_chain*(N_mono+N_cross*N_dangle)) + i*AR+j +2);
        }
    }


    /// Output bead postions to txt file
    int Nang_tot_nw, Nang_perchain;
    out_config.open(fout_conf, ios::out);
    out_config << "\\ Initial Configuration - Tsai-Wei Lin\n" << endl;
    sprintf(null, "# Seed = %ld\n", seed_num);
    out_config << null << endl;
    sprintf(null, "%d atoms", N_chain*(N_mono+N_cross*N_dangle) + Nrod*AR);
    out_config << null << endl;
    sprintf(null, "%d bonds", N_chain*(N_mono-1 + N_cross*N_dangle) + Nrod*(AR-1));
    out_config << null << endl;
    if (N_dangle==0){
        sprintf(null, "%d angles", N_chain*(N_mono-2 + N_cross*(N_dangle))+ Nrod*(((AR-2) < 0) ? 0 : AR-2));
        Nang_tot_nw = N_chain*(N_mono-2 + N_cross*(N_dangle));
        Nang_perchain = N_mono-2 + N_cross*(N_dangle);
    }
    else{
        sprintf(null, "%d angles", N_chain*(N_mono-2 + N_cross*(N_dangle+1))+ Nrod*(((AR-2) < 0) ? 0 : AR-2));  
        Nang_tot_nw = N_chain*(N_mono-2 + N_cross*(N_dangle+1));
        Nang_perchain = N_mono-2 + N_cross*(N_dangle+1);
    }

    out_config << null << endl;
    sprintf(null, "%d dihedrals", 0);
    out_config << null << endl;
    sprintf(null, "%d impropers", 0);
    out_config << null << endl;
    sprintf(null, "%d atom types", 6);
    out_config << null << endl;
    sprintf(null, "%d bond types", 4);
    out_config << null << endl;
    sprintf(null, "%d angle types", 4);
    out_config << null << endl;
    sprintf(null, "%f %f xlo xhi", -box_length/2,box_length/2);
    out_config << null << endl;
    sprintf(null, "%f %f ylo yhi", -box_length/2,box_length/2);
    out_config << null << endl;
    sprintf(null, "%f %f zlo zhi", -box_length/2,box_length/2);
    out_config << null << endl;
    out_config << "\nMasses\n" << endl;
    for (int i=0; i<6; ++i){
        sprintf(null, "%d %f", i+1, 1.0);
        out_config << null << endl;
    }
    
    out_config << "\nAtoms\t# atom-ID mol-ID atom-type q x y z\n" << endl;
    for (int i=0; i<N_chain; ++i){
        for (int j=0; j<chainrx[i].size(); ++j){
            sprintf(null, "%d\t %d\t %d\t %lf\t %lf\t %lf\t %lf", (i*(N_mono+N_cross*N_dangle) + j)+1, i+1, atype[i][j], charge[i][j], chainrx[i][j], chainry[i][j], chainrz[i][j]);
            out_config << null << endl;
        }
    }

    for (int i=0; i<Nrod; ++i){
        for (int j=0; j<AR; ++j){
            sprintf(null, "%d\t %d\t %d\t %lf\t %lf\t %lf\t %lf", N_chain*(N_mono+N_cross*N_dangle) + (i*AR+j+1), N_chain+i+1, patype[i][j], pene_q[i][j], penerx[i][j], penery[i][j], penerz[i][j]);
            out_config << null << endl;
        }
    }

    out_config << "\n\nBonds\n" << endl;
    for (int i=0; i<N_chain; ++i){
        for (int j=0; j<btype[i].size(); ++j){
            sprintf(null, "%d\t %d\t %d\t %d", (i*(N_mono+N_cross*N_dangle) + j)+1-i, btype[i][j], bfirst[i][j] +1, bsecond[i][j] +1);
            out_config << null << endl;
        }
    }

    for (int i=0; i<Nrod; ++i){
        for (int j=0; j<pbtype[i].size(); ++j){
            sprintf(null, "%d\t %d\t %d\t %d", N_chain*(N_mono-1 + N_cross*N_dangle) + (i*AR+j)+1-i, 4, pbfirst[i][j] +1, pbsecond[i][j] +1);
            out_config << null << endl;
        }
    }

    out_config << "\n\nAngles\n" << endl;
    for (int i=0; i<N_chain; ++i){
        for (int j=0; j<angtype[i].size(); ++j){
            sprintf(null, "%d\t %d\t %d\t %d\t %d", (i*(N_mono+N_cross*N_dangle)+j)+1-2*i, angtype[i][j], angfirst[i][j] +1, angsecond[i][j] +1, angthird[i][j] +1);
            out_config << null << endl;
        }
    }

    for (int i=0; i<Nrod; ++i){
        for (int j=0; j<pangtype[i].size(); ++j){
            sprintf(null, "%d\t %d\t %d\t %d\t %d", Nang_tot_nw + (i*AR+j)+1-2*i, pangtype[i][j], pangfirst[i][j] +1, pangsecond[i][j] +1, pangthird[i][j] +1);
            out_config << null << endl;
        }
    }

    free(idum);

    cout << "done!" << endl;
    return 0;
}


//////////////////// FUNCTIONS ///////////////////////////////////////

/*ran1() - Return a random floating point value between 0.0 and
   1.0 exclusive.  If idum is negative, a new series starts (and
   idum is made positive so that subsequent calls using an unchanged
   idum will continue in the same sequence).
   http://cse.unl.edu/~choueiry/Code/urbcsp.c
 */
double ran1(long *idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if(*idum <= 0)
    {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for(j=NTAB+7;j>=0;--j)
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if(*idum<0) *idum+=IM1;
            if(j<NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum<0) *idum += IM1;
    k=idum2/IQ2;
    if(*idum<0) idum2+= IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if(iy<1) iy += IMM1;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

double gasdev(long *idum)
{
    double ran1(long *idum);
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;
    
    if (*idum < 0) iset = 0;
    if (iset == 0)
    {
        do
        {
            v1 = 2.0*ran1(idum)-1.0;
            v2 = 2.0*ran1(idum)-1.0;
            rsq = v1*v1+v2*v2;
        }
        while(rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}

long initRan()
{
   //This will hopefully allow us to have a unique seed even if executed multiple times a second-Got from Mike
   //http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
   unsigned long a = clock();
   unsigned long b = time(NULL);
   unsigned long c = getpid();
   a=a-b;  a=a-c;  a=a^(c >> 13);
   b=b-c;  b=b-a;  b=b^(a << 8);
   c=c-a;  c=c-b;  c=c^(b >> 13);
   a=a-b;  a=a-c;  a=a^(c >> 12);
   b=b-c;  b=b-a;  b=b^(a << 16);
   c=c-a;  c=c-b;  c=c^(b >> 5);
   a=a-b;  a=a-c;  a=a^(c >> 3);
   b=b-c;  b=b-a;  b=b^(a << 10);
   c=c-a;  c=c-b;  c=c^(b >> 15);
   return c%1000000000; //careful here.  Another 0 might break the ran1 (long long instead of just long)
}


