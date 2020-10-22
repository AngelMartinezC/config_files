
/* 
 * Calculate Density, Gamma1 and Sound Speed from thermodynamical
 * pressure and temperature for a 3d data cube (from hhe_thermodynam)
 *
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
//#include <unistd.h>
#include <string.h>
#include "fitsio.h"



void printerror( int status);



int main(){

  char *input1 = "p.fits";
  char *input2 = "T.fits";
  char *output1 = "output1.fits";
  char *output2 = "output2.fits";
  char *output3 = "output3.fits";
  system("rm -rf output*");
  
  
  /* Load fits variables from CFITSIO */
  fitsfile *afptr, *bfptr, *outfptr1, *outfptr2, *outfptr3; 
  int status = 0;    /* CFITSIO status starting from zero */
  int atype, btype, anaxis, bnaxis, check = 1, ii, op;
  long npixels = 1, firstpix[3] = {1,1,1}, ntodo;
  long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1};
  double *apix, *bpix, *cpix, value;
  int image2=1;


  /* Open files and get dimensions */
  fits_open_file(&afptr, input1, READONLY, &status);
  fits_open_file(&bfptr, input2, READONLY, &status);

  fits_get_img_dim(afptr, &anaxis, &status); /* Get dimensions */
  if (image2) fits_get_img_dim(bfptr, &bnaxis, &status);
  fits_get_img_size(afptr, 3, anaxes, &status); /* Real NAXIS */
  if (image2) fits_get_img_size(bfptr, 3, bnaxes, &status);
  
  
  
  /* Calculate thermodynamical density, gamma1 and sound speed from 
   * hhe_thermodynam with total presure and temperature as inputs, and
   * return 3 thermodynamical datacubes */ 
  
  if (check && !fits_create_file(&outfptr1, output1, &status)){
    if (check && !fits_create_file(&outfptr2, output2, &status)){
      if (check && !fits_create_file(&outfptr3, output3, &status)){
      
        fits_copy_header(afptr, outfptr1, &status);
        fits_copy_header(afptr, outfptr2, &status);
        fits_copy_header(afptr, outfptr3, &status);
        
        npixels = anaxes[0];
        
        /* Allocation of data in memory of a row (need to be read as x) */
        apix = (double *) malloc(npixels * sizeof(double)); 
        cpix = (double *) malloc(npixels * sizeof(double));
        if (image2) bpix = (double *) malloc(npixels * sizeof(double));
        
        for (firstpix[2] = 1; firstpix[2] <= anaxes[2]; firstpix[2]++){
          for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++){      
            if (fits_read_pix(afptr, TDOUBLE, firstpix, npixels, NULL, apix,
                              NULL, &status)) break;
            if (image2 && fits_read_pix(bfptr, TDOUBLE, firstpix, npixels,
				          NULL, bpix, NULL, &status)) break;
            
            for(ii=0;ii<npixels; ii++){  /* Iteration over x value (ROW) */ 
              double s1 = apix[ii];      /* Pressure */
              double s2 = bpix[ii];      /* Temperature */
              //printf("%f\n",s1);
              
              /* Create bash command to concatenate with c  */
              char Pre[12], Temp[12];
              sprintf(Pre, "%d ", (int)s1);
              sprintf(Temp, "%d ", (int)s2);
              char command[600];
              strcpy(command," A=`echo ");
              strcat(command,Pre);
              strcat(command,Temp);
              strcat(command," | hhe_thermodynam");
              strcat(command," | egrep 'sound|density|Gamma1'");
              strcat(command," | awk '{print $3}'`; B='_'; echo $A $B");
              strcat(command," | awk '{print $1 $4 $2 $4 $3}'");
              
              //printf("%s    \n",command);
              
              /* Run bash script and read into thermodynamical variables */
              char cadena[100];
              FILE *fp = popen(command, "r"); /* popen to get the bash output */
              fscanf(fp, "%s", cadena);
              char delim[] = "_";
              char *ptr = strtok(cadena, delim);
              char *eptr;
              long double thermo[3]; /* density, Gamma1, sound_speed */
              int count = 0;
              while(ptr != NULL){
		            //printf("'%s'\n", ptr);
		            thermo[count] = strtold(ptr,&eptr);
		            ptr = strtok(NULL, delim);
		            count+=1;
	            }
              pclose(fp);
              
              /* Print "Time step" */
              printf("x: %d    y: %ld    z: %ld\n",ii,firstpix[1],firstpix[2]);
              
              apix[ii] = thermo[0];  /* Density      */
              bpix[ii] = thermo[1];  /* Gamma1       */
              cpix[ii] = thermo[2];  /* Sound speed  */
              }
            fits_write_pix(outfptr1, TDOUBLE, firstpix, npixels, apix, &status);
            fits_write_pix(outfptr2, TDOUBLE, firstpix, npixels, bpix, &status);
            fits_write_pix(outfptr3, TDOUBLE, firstpix, npixels, cpix, &status);
          }
        }
        fits_close_file(outfptr1, &status);
        fits_close_file(outfptr2, &status);
        fits_close_file(outfptr3, &status);
        free(apix);
        free(cpix);
        if (image2) free(bpix);
        printf("%ld  %ld  %ld\n",anaxes[0],anaxes[1],anaxes[2]);
      }
    }
  }

  fits_close_file(afptr, &status);  
  
  if (status) fits_report_error(stderr, status); // print any error message 
  if (image2) fits_close_file(bfptr, &status);
  

  return 0;

}









void printerror( int status)
{
    /* Print cfitsio error messages and exit program */


    if (status){
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
    
}





