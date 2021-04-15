/* 
 * This program calculates the Time Distance relation for a given solar
 * model S of Christensen Dalsgaard. It is performed via Ray-Path 
 * approximation (using ray theory, of course), for HMI's frequency
 * spectral ranges (1e-3,11e-3) mHz. 
 * 
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <iomanip>


#define I     2482          // Rows of file ModelS
#define J     6             // Columns of file ModelS
#define PI    3.1415926535  
#define RSun  69.634e9      // cm


using namespace std;


double Lamb(int l);
double tau(double c, double L, double rad, double w);
double Delta(double c, double L, double rad, double w, double wc);
void Read_variables(void);

double rad[I], c[I], rho[I], p[I], Gamma[I], T[I], wc[I];




// *******************************************************


int main (int argc, char** argv)
{   
    // Initialize variables for modeling
    Read_variables();

    double res_t;
    double res_D;
    int count = 0;

    
    for (int l=200; l<=4001; l+=200){

        ostringstream namestring;
        namestring << l;
        string s(namestring.str());
        string name = "function_l/T_D_" +s+".txt";

        ofstream out (name);

        for (double w=2*PI*11e-3; w>=2*PI*1e-3; w-=2*PI*1e-4){

            double tau_area = 0.0;
            double result_t = 0.0;
            double Del_area = 0.0;
            double result_D = 0.0;

            for (int r=I-2; r>=1; r--){

                res_t = tau(c[r], Lamb(l), rad[r], w);
                res_D = Delta(c[r], Lamb(l), rad[r], w, wc[r]);

                if ((isnan(res_t)==1) and (isnan(res_D)==1)){
                    tau_area += 0.0;
                    Del_area += 0.0;
                }

                else {
                    double b = (rad[r-1]-rad[r])/RSun;
                    if (isnan(result_t)){
                        result_t = 0;
                        result_D = 0;}
                    /*if (res_D/1e8>=60000.0){ // Old bad threshold 
                        res_D = 0;}*/
                    double h_t = res_t + result_t;
                    double h_D = res_D + result_D;
                    double area_t = b*h_t/2;
                    double area_D = b*h_D/2;
                    tau_area += area_t;
                    Del_area += area_D;
                }

                result_t = res_t;
                result_D = res_D;
                
                //out << rad[r] << " " << res_D/1e8 << "\n";
            } 
            out << tau_area/60 << " " << Del_area/1e8 << endl;
        } 
    out.close();

    }

    return 0;

}



// *******************************************************


// Lamb frequency

double Lamb(int l){
    double L = pow(l*(l+1),0.5);
    return L;
}



/* Functions tau and Delta for integration. These functions are calculated
 * via ray-path approximation.
 */

double tau(double c, double L, double rad, double w){
    double res = 2*RSun/(c*pow(1-pow(L*c/(rad*w),2),0.5));
    return res;
}

double Delta(double c, double L, double rad, double w, double wc){
    double res = 2*RSun/(pow(pow(w*rad/(L*c),2) - 1,0.5));
    //res = 2*(c*RSun/rad)/(pow(pow(w/L,2)-pow(wc/L,2)-pow(c/rad,2),0.5));
    return res;
}



/* void function Read_variables. This reads the thermodynamical variables
 * of the model S and put in arrays. It is calculated the density scale
 + height (H) for the cut-off frequency (both in approximation an model).
 */

void Read_variables(void){

    ifstream inFile;
    double C[I][J];

    inFile.open("modelS.txt", ios::in);
    if (! inFile) {
        cerr << "unable to open file C.txt for reading. Exiting" << endl;
        return;}

    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            inFile >> C[i][j];
        }    
    }
    inFile.close();

    for (int i=0; i<I; i++){
            rad[i]   = C[i][0]*RSun;
            c[i]     = C[i][1];
            rho[i]   = C[i][2];
            p[i]     = C[i][3];
            Gamma[i] = C[i][4];
    }

    double H[I];
    double der_log_rho[I];
    double log_rho[I];
    double der_H[I];
       
    for (int i=0; i<I; i++){
        log_rho[i] = log(rho[i]);}

    for (int i=0; i<I; i++){
        if (i == 0){
            der_log_rho[i] = (log_rho[i+1]-log_rho[i])/(rad[i+1]-rad[i]);}
        else if (i == I-1){
            der_log_rho[i] = (log_rho[i]-log_rho[i-1])/(rad[i]-rad[i-1]);}
        else {
            der_log_rho[i] = (log_rho[i+1]-log_rho[i-1])/(2*(rad[i+1]-rad[i]));
        }
        H[i] = -1.0*pow(der_log_rho[i],-1.0);
        wc[i] = c[i]/(2*H[i]);
    }

    /* This computes the cut-off frequency with the given approximation of
     * the derivative of the density scale height. As the derivative of 
     * the natural logarithm of the density is by midpoint mode, there is
     * no guaranty of the derivative of H to work as a smooth function.
     */

    for (int i=0; i<I; i++){
        if (i == 0){
            der_H[i] = (H[i+1]-H[i])/(rad[i+1]-rad[i]);}
        else if (i == I-1){
            der_H[i] = (H[i]-H[i-1])/(rad[i]-rad[i-1]);}
        else {
            der_H[i] = (H[i+1]-H[i-1])/(2*(rad[i+1]-rad[i]));
        }
        //wc[i] = (c[i]/(2*H[i]))*pow(1-2*der_H[i],0.5);
    }

}


