/* Code to make a "HI" word with some selected points around to calculate
 * a postel-projected image of it.
 *
 * The original image is BITPIX 32, which, for a 4096x4096 array, will have 
 * a size around 64 MB (not attached). This is a typical SDO/HMI's image of
 * full-disk magnetic field. The postel projected image projects, with a
 * software, around 3459, 1583 px (LON/LAT Carrington HGC).
 */

#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"

#define PI 3.14159264
#define A 4096

void make_template(short array[A][A]);
void create_template(long array[A][A]);
void open_fits(void);



// ------------------------------------------------------------------
int main(int argc, char *argv[]){

  //int array_new[4][3] = {{0,1,2},{3,4,5},{6,7,8},{9,10,11}};
  //printf("%d\n",array_new[3][2]);

  FILE *data;
  //int x = 900, y = 2400;
  int x = 3459;
  int y = 1583;
  static long array[A][A];


  // -- Make a circle for test
  for (int ii=0; ii<4096; ii++){
    for (int jj=0; jj<4096; jj++){
      int A1 = ii-2040;
      int B1 = jj-2049;
      if (A1*A1+B1*B1>1935*1935){
        // && A1*A1+B1*B1>1650*1650){
        array[jj][ii] = 1.0/0.0;}
      //else if (A1*A1+B1*B1<1500*1500){
        //array[jj][ii] = 0;}
      //else{
        //array[jj][ii] = 1;
    }
  }
 

  // -- Make the word "HI" and save it to a FITS file
  // ---------------------------------------------
  
  for (int ii=0; ii<=35; ii++){
    for (int jj=0; jj<=210; jj++){
      array[jj+y][ii+x-100] = 1;
      array[jj+y][ii+x] = 1;
      array[jj+y][ii+x+150] = 1;
    }
  }

  for (int ii=0; ii<=1; ii++){
    for (int jj=0; jj<=1; jj++){
      array[jj+y-200][ii+x+5] = 1;
      array[jj+y-200][ii+x+10] = 1;
      array[jj+y-200][ii+x+20] = 1;
      array[jj+y-200][ii+x+110] = 1;
      array[jj+y-200][ii+x+22] = 1;
      array[jj+y-200][ii+x+25] = 1;
      array[jj+y-200][ii+x+28] = 1;
      array[jj+y-200][ii+x+31] = 1;
      array[jj+y-200][ii+x+33] = 1;
      array[jj+y-203][ii+x+25] = 1;
      array[jj+y-197][ii+x+25] = 1;
    }
  }

  for (int ii=0; ii<=5; ii++){
    for (int jj=0; jj<=5; jj++){
      array[jj+y-100][ii+x] = 1;
      array[jj+y-100][ii+x+100] = 1;
      array[jj+y-80][ii+x-80] = 1;
      array[jj+y-50][ii+x+110] = 1;
      array[jj+y-80][ii+x+110] = 1;
    }
  }

  for (int ii=0; ii<=20; ii++){
    for (int jj=0; jj<=20; jj++){
      array[jj+y-50][ii+x-50] = 1;
      array[jj+y-90][ii+x-50] = 1;
      array[jj+y-50][ii+x-90] = 1;
      array[jj+y-50][ii+x+90] = 1;
      array[jj+y-50][ii+x+40] = 1;
    }
  }

  for (int ii=0; ii<=20; ii++){
    for (int jj=0; jj<=20; jj++){
      array[jj+y-50][ii+x-50] = 1;
      array[jj+y-90][ii+x-50] = 1;
      array[jj+y-50][ii+x-90] = 1;
      array[jj+y-50][ii+x+90] = 1;
      array[jj+y-50][ii+x+40] = 1;
    }
  }

  for (int i=0; i<=100; i++){
    for (int j=0; j<=35; j++){
      array[j+y+125-50+17][i+x-100] = 1; //Horizontal H
      array[j+y][i+x+150-50+17] = 1;
      array[j+y+210-52+17][i+x+150-50+17] = 1;
    }
  }
  
  // -- Create a point in center (reference point)
  
  for (int i=0; i<=100; i++){
    for (int j=0; j<=10; j++){
      array[j+2024+10][i+2024-50] = 1;//4547.0490723;
      array[j+2024][i+2024-50] = 1;
      array[j+2024-10][i+2024-50] = 1;//-4547.0490723;
    }
  }
  
  // -- Saving to fits file from cfitsio program
  create_template(array);
  
  //open_fits(); 

}

// ------------------------------------------------------------------


/* Make a simple file
 */
void make_template(int array[A][A]){
  
  FILE *data;

  // -- Create a TEMPLATE2
  data = fopen("TEMPLATE2.txt","w");
  for (int i=0; i<A; i++){
    for (int j=0; j<A; j++){
      if (j==A-1)
        fprintf(data, "%d\n",array[j][i]);
      else  
        fprintf(data, "%d ",array[j][i]);
    }
  }
  fclose(data);
}



/* Function to make fits files from array data. By default, as for the
 * further requirements in working with long data format, I set the datatype
 * to be LONG_IMG (BITPIX 32). 
 */
void create_template(long array[A][A]){
 
  // -- pointer to the FITS file; defined in fitsio.h 
  fitsfile *fptr;        

  int status, ii, jj;
  double fpixel = 1, naxis = 2, nelements;
  long naxes[2] = { A, A};  

  system("rm TEMPLATE2.fits");
  status = 0;         // -- initialize status before calling routines 
  fits_create_file(&fptr, "TEMPLATE2.fits", &status);   
  
  // -- Write a keyword; must pass the ADDRESS of the value 
  fits_create_img(fptr, LONG_IMG, naxis, naxes, &status);
  nelements = naxes[0] * naxes[1];   // -- number of pixels to write
  
  // -- Write the array of integers to the image 
  fits_write_img(fptr, TLONG, fpixel, nelements, array[0], &status);
  fits_close_file(fptr, &status);            // -- close the file 
  fits_report_error(stderr, status);  // -- print out any errors 

}


/* Checkout 
 */
void open_fits(void){
  int status = 0;
  fitsfile *fptr;
  fits_open_file(&fptr, "SAMPFILE.fits", READWRITE, &status);

}
