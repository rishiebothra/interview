
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<dirent.h>
#include"jsmn.h"

/****************************Global Declarations************************/
#define JSON_FILE_PATH "./input_data/lbmPara.json"
#define BUFFER_SIZE 5000
#define MAX_TOKEN_COUNT 128
#define QMAX    19
#define exLength  5

int exLength2 = exLength*2;
int IMAX,JMAX,KMAX,TMAX;
double tau,nu,force,fextx,fexty,fextz,pixSize,tolerance;
char PBC_switch[128];
double rho_i,rho_o;
char flow_direction[128];
char periodic_x[128],periodic_y[128],periodic_z[128];
char vel_output[128];
char collision_model[128];
char input_file_name[128],output_file_name[128];
char dataFormat[128];
char tmpStr[128];
char inName[128] = "./";
char outName[128] = "./";

double ****f,****f_coll; //distribution function before and after collision
double ***ux,***uy,***uz; //macroscopic velocity
double ***ux_old,***uy_old,***uz_old; //velocity at previous ts
double ***rho,***press; //density and pressure
double ***rho_old,***press_old; //density and pressure at previous ts
int ***isnode; //isnode=0: pore; isnode=1: solid
double ***Fx,***Fy,***Fz; //external forces
int    ex[QMAX],ey[QMAX],ez[QMAX]; //particle velocities
double eex[QMAX],eey[QMAX],eez[QMAX]; //particle velocities
int    opp[QMAX]; //bounce-back direction
double w[QMAX]; //discrete weighting factors
double SCAL[QMAX]; //Scaling factors //MRT
double omg[QMAX],fv[QMAX],m[QMAX],Sm[QMAX],cm[QMAX];
double cv[QMAX]; //f change due to collision
double cs,cs2,acs2,acs4,acs4_2,acs2_2; //speed of sound and its variants
int    ntime; //time steps
double error,error_old,error_diff;
clock_t tStart,tEnd,tUsed;

double Gamhat[QMAX],Lambdahat[QMAX];

FILE   *fpo_vel; //output velocity fields
FILE   *fpo_sum; //output summary
FILE   *fpo_K;   //output Kx
/****************************Global Declarations************************/

void Initialize_f(void);     //DONE
void Initialize_macro(void); //DONE
void kernel(void);           //DONE //Collision+Streaming+Macro
void stream(void);           //DONE
void macro(void);            //DONE
void convergence(void);      //DONE
void results(void);          //DONE
void free_memory(void);      //DONE
void collision_SRT(void);    //DONE
void collision_MRT(void);    //DONE
void PBCondition(void);      //XYZ directions
void lbmPara(void);          //DONE
void geoDef(void);           //DONE
void cvtInt(char *str,char *inP);
double Meq(double RHO,double JX,double JY,double JZ,double JDOTJ,int nq);
double FtoM(double f0,double f1,double f2,double f3,double f4,double f5,double f6,double f7,double f8,double f9,double f10,double f11,double f12,double f13,double f14,double f15,double f16,double f17,double f18,int nq);
double MtoF(double cm0,double cm1,double cm2,double cm3,double cm4,double cm5,double cm6,double cm7,double cm8,double cm9,double cm10,double cm11,double cm12,double cm13,double cm14,double cm15,double cm16,double cm17,double cm18,int nq);

void readfile(char* filepath, char* fileContent);
int parseJSON(char *filepath, void callback(char *, char*));
void mycallback(char *key, char* value);
void LBM_calculation(void);

int main(int argc, char* argv[]){


  strcat(inName,argv[1]);
  strcat(inName,"/input-image.dat");


  if(strcmp(argv[1],"input_data")!=0)
    {
      printf(" input folder name!\n");

      //return 0;
    }
  else if(strcmp(argv[2],"output_data")!=0){
    printf(" output folder name!\n");
    //return 0;
  }
  else{



    tStart = clock();

    //Reading LBM Input Parameters from json file
    parseJSON(JSON_FILE_PATH, mycallback);

    int dirU;
    if(strcmp(flow_direction,"all")==0){
      for(dirU=0;dirU<=2;dirU++){

      if(dirU==0){

	strcpy(flow_direction,"x");
	strcpy(periodic_x,"on");
	IMAX += exLength2;
	JMAX += 2;
	KMAX += 2;
	fextx = force;
      }

      if(dirU==1){

	strcpy(flow_direction,"y");
	strcpy(periodic_y,"on");
	JMAX += exLength2;
	IMAX += 2;
	KMAX += 2;
	fexty = force;
      }

      if(dirU==2){

	strcpy(flow_direction,"z");
	strcpy(periodic_z,"on");
	KMAX += exLength2;
	JMAX += 2;
	IMAX += 2;
	fextz = force;
      }



      strcpy(tmpStr, flow_direction);
      strcat(tmpStr, output_file_name);
      strcpy(output_file_name, tmpStr);
      strcat(outName,argv[2]);
      strcat(outName,"/");
      strcat(outName,output_file_name);

      LBM_calculation();

     //Reset the input parameters
      IMAX = JMAX = KMAX = 100;
      fextx = fexty = fextz = 0.0;
      strcpy(periodic_x,"off");
      strcpy(periodic_y,"off");
      strcpy(periodic_z,"off");

      //reset output_file_name to "Perm_LBM.json"
      strcpy(output_file_name,"Perm_LBM.json");
      strcpy(outName,"./");

      } //for loop dirU

    } //if flow_direction == all
    else if(strcmp(flow_direction,"x")==0 || strcmp(flow_direction,"y")==0 || strcmp(flow_direction,"z")==0){

      if(strcmp(flow_direction,"x")==0){IMAX += exLength2;JMAX += 2;KMAX += 2;fextx = force;}
      if(strcmp(flow_direction,"y")==0){JMAX += exLength2;IMAX += 2;KMAX += 2;fexty = force;}
      if(strcmp(flow_direction,"z")==0){KMAX += exLength2;JMAX += 2;IMAX += 2;fextz = force;}


      strcpy(tmpStr, flow_direction);
      strcat(tmpStr, output_file_name);
      strcpy(output_file_name, tmpStr);

      strcat(outName,argv[2]);
      strcat(outName,"/");
      strcat(outName,output_file_name);

      LBM_calculation();

    }
  else printf("Flow Direction Can only be x, y, z, or all!\n");

  }//if input folder argument not match!

  return 0;

}

void LBM_calculation(void){

  printf("----------------------------------------------------------------\n");
  printf("Flow Direction :  %s\n",flow_direction);
  printf("IMAX =      %d\n",IMAX);
  printf("JMAX =      %d\n",JMAX);
  printf("KMAX =      %d\n",KMAX);
  printf("tau  =      %.4f\n",tau);
  printf("nu   =      %.4f\n",nu);
  printf("external_force_x        =  %.4e\n",fextx);
  printf("external_force_y        =  %.4e\n",fexty);
  printf("external_force_z        =  %.4e\n",fextz);
  printf("convergence_tolerance   =  %.4e\n",tolerance);
  printf("pixel_resolution(meter) =  %.4e\n",pixSize);
  printf("periodic_x :      %s\n",periodic_x);
  printf("periodic_y :      %s\n",periodic_y);
  printf("periodic_z :      %s\n",periodic_z);
  printf("pressure_BC:      %s\n",PBC_switch);
  printf("inlet_rho:        %f\n",rho_i);
  printf("outlet_rho:       %f\n",rho_o);
  printf("vel_output :      %s\n",vel_output);
  printf("collision_model:  %s\n",collision_model);
  printf("dataFormat :      %s\n",dataFormat);
  printf("input_file_name:  %s\n",input_file_name);
  printf("output_file_name: %s\n",output_file_name);
  printf("----------------------------------------------------------------\n");


  //Define LBM intrinsic parameters
  lbmPara();


  //Read and define solid boundaries
  geoDef();

  // Initilization of macroscopic and distribution quantities
  ntime = 0;

  Initialize_macro();
  Initialize_f();

  // Time-involving
  for(ntime=1;ntime<=TMAX;ntime++)
    {

      kernel();


      if(ntime%10==0){

	printf("Time Step = %d  ",ntime);

	convergence();

	if((error<=tolerance) || (ntime==TMAX)){

	  printf("Solution is converged sucessifully!\n");

	  results();

	  break;
	}
      }

    }


  free_memory();


  if(strcmp(vel_output,"on") == 0)  fclose(fpo_vel);


  tEnd = clock();
  double tUsed = tEnd - tStart;
  printf ("Time used: (%f seconds).\n",((float)tUsed)/CLOCKS_PER_SEC);


}


void lbmPara(void){

  // Set particle velocities D3Q19
  eex[0] =0.0;  eey[0] =0.0;	eez[0] =0.0;
  eex[1] =+1.0; eey[1] =0.0;	eez[1] =0.0;
  eex[2] =-1.0; eey[2] =0.0;	eez[2] =0.0;
  eex[3] =0.0;  eey[3] =+1.0;	eez[3] =0.0;
  eex[4] =0.0;  eey[4] =-1.0;	eez[4] =0.0;
  eex[5] =0.0;  eey[5] =0.0;	eez[5] =+1.0;
  eex[6] =0.0;  eey[6] =0.0;	eez[6] =-1.0;
  eex[7] =+1.0; eey[7] =+1.0;	eez[7] =0.0;
  //  eex[8] =-1.0; eey[8] =-1.0;	eez[8] =0.0;
  eex[8] =-1.0; eey[8] =+1.0;   eez[8] = 0.0;
  eex[9] =+1.0; eey[9] =-1.0;	eez[9] =0.0;
  //  eex[10] =-1.0;eey[10] =+1.0;	eez[10] =0.0;
  eex[10] =-1.0;eey[10] =-1.0;  eez[10] = 0.0;
  eex[11] =+1.0;eey[11] =0.0;	eez[11] =+1.0;
  //  eex[12] =-1.0;eey[12] =0.0;	eez[12] =-1.0;
  eex[12] =-1.0;eey[12] =0.0;	eez[12] =+1.0;
  //  eex[13] =-1.0;eey[13] =0.0;	eez[13] =+1.0;
  eex[13] =+1.0;eey[13] =0.0;	eez[13] =-1.0;
  //eex[14] =+1.0;eey[14] =0.0;	eez[14] =-1.0;
  eex[14] =-1.0;eey[14] =0.0;	eez[14] =-1.0;
  eex[15] =0.0; eey[15] =+1.0;	eez[15] =+1.0;
  //  eex[16] =0.0; eey[16] =-1.0;	eez[16] =-1.0;
  eex[16] =0.0; eey[16] =-1.0;	eez[16] =+1.0;
  eex[17] =0.0; eey[17] =+1.0;	eez[17] =-1.0;
  //  eex[18] =0.0; eey[18] =-1.0;	eez[18] =+1.0;
  eex[18] =0.0; eey[18] =-1.0;	eez[18] =-1.0;


  ex[0] =0;  ey[0] =0;	ez[0] =0;
  ex[1] =+1; ey[1] =0;	ez[1] =0;
  ex[2] =-1; ey[2] =0;	ez[2] =0;
  ex[3] =0;  ey[3] =+1;	ez[3] =0;
  ex[4] =0;  ey[4] =-1;	ez[4] =0;
  ex[5] =0;  ey[5] =0;	ez[5] =+1;
  ex[6] =0;  ey[6] =0;	ez[6] =-1;
  ex[7] =+1; ey[7] =+1;	ez[7] =0;
  ex[8] =-1; ey[8] =+1; ez[8] = 0;
  ex[9] =+1; ey[9] =-1;	ez[9] =0;
  ex[10] =-1;ey[10] =-1;ez[10] = 0;
  ex[11] =+1;ey[11] =0;	ez[11] =+1;
  ex[12] =-1;ey[12] =0;	ez[12] =+1;
  ex[13] =+1;ey[13] =0;	ez[13] =-1;
  ex[14] =-1;ey[14] =0;	ez[14] =-1;
  ex[15] =0; ey[15] =+1;ez[15] =+1;
  ex[16] =0; ey[16] =-1;ez[16] =+1;
  ex[17] =0; ey[17] =+1;ez[17] =-1;
  ex[18] =0; ey[18] =-1;ez[18] =-1;



  w[0]  = 1.0/3.0;
  w[1]  = 1.0/18.0;
  w[2]  = 1.0/18.0;
  w[3]  = 1.0/18.0;
  w[4]  = 1.0/18.0;
  w[5]  = 1.0/18.0;
  w[6]  = 1.0/18.0;
  w[7]  = 1.0/36.0;
  w[8]  = 1.0/36.0;
  w[9]  = 1.0/36.0;
  w[10] = 1.0/36.0;
  w[11] = 1.0/36.0;
  w[12] = 1.0/36.0;
  w[13] = 1.0/36.0;
  w[14] = 1.0/36.0;
  w[15] = 1.0/36.0;
  w[16] = 1.0/36.0;
  w[17] = 1.0/36.0;
  w[18] = 1.0/36.0;


  opp[1] = 2;
  opp[2] = 1;
  opp[3] = 4;
  opp[4] = 3;
  opp[5] = 6;
  opp[6] = 5;
  opp[7] = 10;
  opp[8] = 9;
  opp[9] = 8;
  opp[10]= 7;
  opp[11]= 14;
  opp[12]= 13;
  opp[13]= 12;
  opp[14]= 11;
  opp[15]= 18;
  opp[16]= 17;
  opp[17]= 16;
  opp[18]= 15;


  cs = sqrt(1.0/3.0);
  cs2 = cs*cs;
  acs2 = 1.0/cs2;
  acs2_2 = 0.5*acs2;
  acs4	 = acs2*acs2;
  acs4_2 = 0.5*acs4;


 // Normalization/scaling factors

  SCAL[0]  = 1.0/19.0;
  SCAL[1]  = 1.0/2394.0;
  SCAL[2]  = 1.0/252.0;
  SCAL[3]  = 1.0/10.0;
  SCAL[4]  = 1.0/40.0;
  SCAL[5]  = 1.0/10.0;
  SCAL[6]  = 1.0/40.0;
  SCAL[7]  = 1.0/10.0;
  SCAL[8]  = 1.0/40.0;
  SCAL[9]  = 1.0/36.0;
  SCAL[10] = 1.0/72.0;
  SCAL[11] = 1.0/12.0;
  SCAL[12] = 1.0/24.0;
  SCAL[13] = 1.0/4.0;
  SCAL[14] = SCAL[13];
  SCAL[15] = SCAL[13];
  SCAL[16] = 1.0/8.0;
  SCAL[17] = SCAL[16];
  SCAL[18] = SCAL[16];


  // Relaxation rates in moment space: omg[nv] = 1/relaxation_time
  omg[0]  = 1.0;      //density             //conserved
  omg[1]  = 1.19;     //internal energy
  omg[2]  = 1.4;      //energy squared
  omg[3]  = 1.0;      //momentum flux in x  //conserved
  omg[4]  = 1.2;      //"energy" flux in x
  omg[5]  = 1.0;      //momentum flux in y  //conserved
  omg[6]  = 1.4;      //"energy" flux in y
  omg[7]  = 1.0;      //momentum flux in z  //conserved
  omg[8]  = 1.8;      //"energy" flux in z
  omg[9]  = 1.0/tau;  //Pxx: symmetric traceless viscous stress tensor
  omg[10] = 1.4;      //Additional moments: second rank tensor
  omg[11] = omg[9];   //Pww=Pyy-Pzz
  omg[12] = 1.4;      //Additional moments: second rank tensor
  omg[13] = omg[9];   //Pxy
  omg[14] = omg[9];   //Pyz
  omg[15] = omg[9];   //Pxz
  omg[16] = 1.98;     //Additional moments: third rank tensor
  omg[17] = 1.98;     //Additional moments: third rank tensor
  omg[18] = 1.98;     //Additional moments: third rank tensor


  //Relaxation times with scaling in moment space
  for(int nv=0;nv<QMAX;nv++)
    {

      Lambdahat[nv]= omg[nv];
      Gamhat[nv]= SCAL[nv]*Lambdahat[nv];
    }

}


void geoDef(void){

  int *one;
  int MAX_LINE;

  if(strcmp(flow_direction,"x")==0) MAX_LINE = (IMAX-exLength2)*(JMAX-2)*(KMAX-2);
  if(strcmp(flow_direction,"y")==0) MAX_LINE = (IMAX-2)*(JMAX-exLength2)*(KMAX-2);
  if(strcmp(flow_direction,"z")==0) MAX_LINE = (IMAX-2)*(JMAX-2)*(KMAX-exLength2);

  /*
  struct dirent *de; //struct dirent *de;  // Pointer for directory entry

  DIR *dr = opendir(".");   // opendir() returns a pointer of DIR type.
  */

  FILE *geo_data;


  isnode = (int ***) malloc(IMAX*sizeof(int **));

  for(int i=0;i<IMAX;i++){
    isnode[i] = (int **) malloc(JMAX*sizeof(int *));

    for(int j=0;j<JMAX;j++){
      isnode[i][j] = (int *)malloc(KMAX*sizeof(int));
    }
  }


  one = (int*)malloc(sizeof(int)*MAX_LINE);

    geo_data = fopen(inName,"r");

    if(strcmp(dataFormat,"DoE")==0)
      {

	if(geo_data == NULL)
	  {
	    printf("Geometry file failed to open!\n");
	  }
	else
	  {
	    printf("Geometry file opened sucessfully!\n");
	    //printf("File Name: %s\n", de->d_name);
	    printf("File Name: %s\n", input_file_name);
	    printf("----------------------------------------------------------------\n");

	    char line[MAX_LINE];
	    int x,y,z,i,j,k;
	    int isd;
	    int counter = 0;

	    do{
	      fscanf(geo_data, "%d%d%d%d", &x,&y,&z,&isd);
	      i = x;
	      j = y;
	      k = z;
	      isnode[i][j][k] = isd;
	      counter++;
	    } while(fgets(line, sizeof(line), geo_data) != NULL);

	    printf("Number of lines in the file is %d\n", counter);
	  }
      }

    else if(strcmp(dataFormat,"Binary_only")==0)
      {

	if(geo_data == NULL)
	  {		printf("Geometry file failed to open!\n");
	  }
	else
	  {
	    printf("Geometry file opened sucessfully!\n");
	    //printf("File Name: %s\n", de->d_name);
	    printf("File Name: %s\n", input_file_name);
	    printf("----------------------------------------------------------------\n");

	    for(int i=0;i<MAX_LINE;i++){

	      fscanf(geo_data,"%i",&one[i]);
	    }

	    if(strcmp(flow_direction,"x")==0){
	      for(int i=exLength;i<IMAX-exLength;i++)
		for(int j=1;j<JMAX-1;j++)
		  for(int k=1;k<KMAX-1;k++)
		    isnode[i][j][k] = one[(JMAX-2)*(KMAX-2)*(i-exLength)+(KMAX-2)*(j-1)+(k-1)];
	    }
	    if(strcmp(flow_direction,"y")==0){
	      for(int i=1;i<IMAX-1;i++)
		for(int j=exLength;j<JMAX-exLength;j++)
		  for(int k=1;k<KMAX-1;k++)
		    isnode[i][j][k] = one[(JMAX-exLength2)*(KMAX-2)*(i-1)+(KMAX-2)*(j-exLength)+(k-1)];

		}
	    if(strcmp(flow_direction,"z")==0){
	      for(int i=1;i<IMAX-1;i++)
		for(int j=1;j<JMAX-1;j++)
		  for(int k=exLength;k<KMAX-exLength;k++)
		    isnode[i][j][k] = one[(JMAX-2)*(KMAX-exLength2)*(i-1)+(KMAX-exLength2)*(j-1)+(k-exLength)];

	    }


	  }
      }
    else printf("Please Enter Correct dataFormat!\n");
    //      }

    //  closedir(dr);

    if(strcmp(flow_direction,"x")==0){
      for(int i=0;i<IMAX;i++)
	for(int j=0;j<JMAX;j++)
	  for(int k=0;k<KMAX;k++){
	    if(i<exLength || i>=IMAX-exLength)   isnode[i][j][k] = 0;
	    //	  isnode[i][j][k] = 0;
	    if(j==0 || j==JMAX-1) isnode[i][j][k] = 1;
	    if(k==0 || k==KMAX-1) isnode[i][j][k] = 1;
	  }
    }

    if(strcmp(flow_direction,"y")==0){
      for(int i=0;i<IMAX;i++)
	for(int j=0;j<JMAX;j++)
	  for(int k=0;k<KMAX;k++){
	    if(j<exLength || j>=JMAX-exLength)   isnode[i][j][k] = 0;
	    //	  isnode[i][j][k] = 0;
	    if(i==0 || i==IMAX-1) isnode[i][j][k] = 1;
	  if(k==0 || k==KMAX-1) isnode[i][j][k] = 1;
	  }

    }

    if(strcmp(flow_direction,"z")==0){
      for(int i=0;i<IMAX;i++)
	for(int j=0;j<JMAX;j++)
	  for(int k=0;k<KMAX;k++){
	    if(k<exLength || k>=KMAX-exLength)   isnode[i][j][k] = 0;
	    //	  isnode[i][j][k] = 0;
	    if(i==0 || i==IMAX-1) isnode[i][j][k] = 1;
	    if(j==0 || j==JMAX-1) isnode[i][j][k] = 1;
	  }

    }


}


void Initialize_macro(void){

  ux        = (double ***) malloc(IMAX*sizeof(double **));
  uy        = (double ***) malloc(IMAX*sizeof(double **));
  uz        = (double ***) malloc(IMAX*sizeof(double **));
  rho       = (double ***) malloc(IMAX*sizeof(double **));
  ux_old    = (double ***) malloc(IMAX*sizeof(double **));
  uy_old    = (double ***) malloc(IMAX*sizeof(double **));
  uz_old    = (double ***) malloc(IMAX*sizeof(double **));
  rho_old   = (double ***) malloc(IMAX*sizeof(double **));
  press     = (double ***) malloc(IMAX*sizeof(double **));
  Fx        = (double ***) malloc(IMAX*sizeof(double **));
  Fy        = (double ***) malloc(IMAX*sizeof(double **));
  Fz        = (double ***) malloc(IMAX*sizeof(double **));

  for(int i=0;i<IMAX;i++){
    ux[i]        = (double **) malloc(JMAX*sizeof(double *));
    uy[i]        = (double **) malloc(JMAX*sizeof(double *));
    uz[i]        = (double **) malloc(JMAX*sizeof(double *));
    rho[i]       = (double **) malloc(JMAX*sizeof(double *));
    ux_old[i]    = (double **) malloc(JMAX*sizeof(double *));
    uy_old[i]    = (double **) malloc(JMAX*sizeof(double *));
    uz_old[i]    = (double **) malloc(JMAX*sizeof(double *));
    rho_old[i]   = (double **) malloc(JMAX*sizeof(double *));
    press[i]     = (double **) malloc(JMAX*sizeof(double *));
    Fx[i]        = (double **) malloc(JMAX*sizeof(double *));
    Fy[i]        = (double **) malloc(JMAX*sizeof(double *));
    Fz[i]        = (double **) malloc(JMAX*sizeof(double *));
  }

  for(int i=0;i<IMAX;i++){
    for(int j=0;j<JMAX;j++){
      ux[i][j]        = (double *)malloc(KMAX*sizeof(double));
      uy[i][j]        = (double *)malloc(KMAX*sizeof(double));
      uz[i][j]        = (double *)malloc(KMAX*sizeof(double));
      rho[i][j]       = (double *)malloc(KMAX*sizeof(double));
      ux_old[i][j]    = (double *)malloc(KMAX*sizeof(double));
      uy_old[i][j]    = (double *)malloc(KMAX*sizeof(double));
      uz_old[i][j]    = (double *)malloc(KMAX*sizeof(double));
      rho_old[i][j]   = (double *)malloc(KMAX*sizeof(double));
      press[i][j]     = (double *)malloc(KMAX*sizeof(double));
      Fx[i][j]        = (double *)malloc(KMAX*sizeof(double));
      Fy[i][j]        = (double *)malloc(KMAX*sizeof(double));
      Fz[i][j]        = (double *)malloc(KMAX*sizeof(double));
    }
  }
  for(int i=0;i<IMAX;i++){
    for(int j=0;j<JMAX;j++){
      for(int k=0;k<KMAX;k++)
	{
	  if(isnode[i][j][k]==0){
	    ux[i][j][k]   = 0.0;
	    uy[i][j][k]   = 0.0;
	    uz[i][j][k]   = 0.0;

	    rho[i][j][k]  = 1.0;
	    press[i][j][k] = rho[i][j][k]*cs2;

	    Fx[i][j][k] = fextx;
	    Fy[i][j][k] = fexty;
	    Fz[i][j][k] = fextz;

	    ux_old[i][j][k]   = ux[i][j][k];
	    uy_old[i][j][k]   = uy[i][j][k];
	    uz_old[i][j][k]   = uz[i][j][k];
	    rho_old[i][j][k]  = rho[i][j][k];

	  }
	}
    }
  }

}

void Initialize_f(void){

  double feq;
  double u_x,u_y,u_z,uu,e_x,e_y,e_z,eu;
  int nv;

  f = (double ****)malloc(IMAX*sizeof(double ***));
  f_coll = (double ****)malloc(IMAX*sizeof(double ***));

  for(int i=0;i<IMAX;i++){
    f[i] = (double ***)malloc(JMAX*sizeof(double **));
    f_coll[i] = (double ***)malloc(JMAX*sizeof(double **));
  }

  for(int i=0;i<IMAX;i++){
    for(int j=0;j<JMAX;j++){
      f[i][j] = (double **)malloc(KMAX*sizeof(double *));
      f_coll[i][j] = (double **)malloc(KMAX*sizeof(double *));
    }
  }

  for(int i=0;i<IMAX;i++){
    for(int j=0;j<JMAX;j++){
      for(int k=0;k<KMAX;k++)
	{
	  f[i][j][k]= (double *)malloc(QMAX*sizeof(double));
	  f_coll[i][j][k]= (double *)malloc(QMAX*sizeof(double));

	  if(isnode[i][j][k]==0){

	    u_x = ux[i][j][k];
	    u_y = uy[i][j][k];
	    u_z = uz[i][j][k];
	    uu = u_x*u_x+u_y*u_y+u_z*u_z;

	    for(int nv=0;nv<QMAX;nv++){

	      e_x = eex[nv];
	      e_y = eey[nv];
	      e_z = eez[nv];
	      eu = e_x*u_x+e_y*u_y+e_z*u_z;

	      feq = w[nv]*rho[i][j][k]*(1.0+acs2*eu+0.5*acs2*acs2*eu*eu-0.5*acs2*uu);
	      f[i][j][k][nv] = feq;
	    }
	  }
	}
    }
  }


}

void kernel(void){

  //  clock_t t0,t1,t2,t3;

  //  t0 = clock();

  if(strcmp(collision_model,"SRT")==0) collision_SRT();
  if(strcmp(collision_model,"MRT")==0) collision_MRT();

  //  t1 = clock();

  stream();
  //  t2 = clock();

  macro();


}

void collision_SRT(void){

  double u_x,u_y,u_z,uu;
  double e_x,e_y,e_z,eu;
  double F_x,F_y,F_z;
  double feq;
  double e_u_dot_F;
  double atau, taufac;
  int i,j,k,is,js,ks;

  atau = 1.0/tau;
  taufac = 1.0 - 0.5*atau;

  for(i=0;i<IMAX;i++)
    for(j=0;j<JMAX;j++)
      for(k=0;k<KMAX;k++)
	{
	  if(isnode[i][j][k]==0){

	    Fx[i][j][k] = fextx;
	    Fy[i][j][k] = fexty;
	    Fz[i][j][k] = fextz;

	    F_x = Fx[i][j][k];
	    F_y = Fy[i][j][k];
	    F_z = Fz[i][j][k];
	    u_x = ux[i][j][k];
	    u_y = uy[i][j][k];
	    u_z = uz[i][j][k];
	    uu = u_x*u_x + u_y*u_y + u_z*u_z;

	    for(int nv=0;nv<QMAX;nv++){

	      e_x = eex[nv];
	      e_y = eey[nv];
	      e_z = eez[nv];
	      eu = e_x*u_x+e_y*u_y+e_z*u_z;
	      e_u_dot_F = (e_x-u_x)*F_x + (e_y-u_y)*F_y + (e_z-u_z)*F_z;

	      feq = w[nv]*rho[i][j][k]*(1.0+acs2*eu+0.5*acs2*acs2*eu*eu-0.5*acs2*uu);
	      cv[nv] = -atau*(f[i][j][k][nv]-feq)+taufac*w[nv]*(acs2*e_u_dot_F+acs2*acs2*e_u_dot_F*eu);

	      f_coll[i][j][k][nv] = f[i][j][k][nv] + cv[nv];

	    }
	  }
	}

}


void stream(void){
  int i,j,k,nv;
  int is,js,ks;


  if(strcmp(PBC_switch,"on")==0)
    {
      for(i=1;i<IMAX-1;i++)
	for(j=1;j<JMAX-1;j++)
	  for(k=1;k<KMAX-1;k++)
	    if(isnode[i][j][k]==0){
	      for(nv=0;nv<QMAX;nv++){
		is = i - ex[nv];
		js = j - ey[nv];
		ks = k - ez[nv];

		// Free-Streaming from source nodes
		if((isnode[i][j][k]==0)&&(isnode[is][js][ks]==0))
		  f[i][j][k][nv] = f_coll[is][js][ks][nv];

		// Boundary Condition for near-solid nodes: Bounce back
		if((isnode[i][j][k]==0)&&(isnode[is][js][ks]==1))
		  f[i][j][k][nv] = f_coll[i][j][k][opp[nv]];
	      }
	    }
      PBCondition();
    }

  else{
    for(i=0;i<IMAX;i++)
      for(j=0;j<JMAX;j++)
	for(k=0;k<KMAX;k++)
	  if(isnode[i][j][k]==0){
	    for(nv=0;nv<QMAX;nv++){
	      is = i - ex[nv];
	      js = j - ey[nv];
	      ks = k - ez[nv];

	      if(strcmp(periodic_x,"on")==0){
		if(is<0)      is = is + IMAX;
		if(is>IMAX-1) is = is - IMAX;
	      }
	      if(strcmp(periodic_y,"on")==0){
		if(js<0)      js = js + JMAX;
		if(js>JMAX-1) js = js - JMAX;
	      }

	      if(strcmp(periodic_z,"on")==0){
		if(ks<0)      ks = ks + KMAX;
		if(ks>KMAX-1) ks = ks - KMAX;
	      }

	      // Free-Streaming from source nodes
	      if((isnode[i][j][k]==0)&&(isnode[is][js][ks]==0))
		f[i][j][k][nv] = f_coll[is][js][ks][nv];

	      // Boundary Condition for near-solid nodes: Bounce back
	      if((isnode[i][j][k]==0)&&(isnode[is][js][ks]==1))
		f[i][j][k][nv] = f_coll[i][j][k][opp[nv]];
	    }

	  }

  }

}

void PBCondition(void){

  double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18;
  double Nx_y,Nx_z,Ny_x,Ny_z,Nz_x,Nz_y;
  double rho_tmp,ux_tmp,uy_tmp,uz_tmp;

  for(int i=0;i<IMAX;i++)
    for(int j=0;j<JMAX;j++)
      for(int k=0;k<KMAX;k++)
	{
	  if(isnode[i][j][k]==0){

	    f0=f[i][j][k][0];
	    f1=f[i][j][k][1];  f2 =f[i][j][k][2];
	    f3=f[i][j][k][3];  f4 =f[i][j][k][4];
	    f5=f[i][j][k][5];  f6 =f[i][j][k][6];
	    f7=f[i][j][k][7];  f8 =f[i][j][k][8]; f9 =f[i][j][k][9]; f10=f[i][j][k][10];
	    f11=f[i][j][k][11];f12=f[i][j][k][12];f13=f[i][j][k][13];f14=f[i][j][k][14];
	    f15=f[i][j][k][15];f16=f[i][j][k][16];f17=f[i][j][k][17];f18=f[i][j][k][18];



	    if(strcmp(flow_direction,"x")==0){
	      //X-Direction_inlet
	      if(i==0){
		rho_tmp = rho_i;
		uy_tmp  = 0.0; //not neccessarily zero
		uz_tmp  = 0.0; //not neccessarily zero
		ux_tmp  = 1.0-1.0/rho_tmp*(f0+f3+f4+f5+f6+f15+f17+f16+f18+2.0*(f2+f8+f10+f12+f14));
		Nx_y = 0.5*(f3+f15+f17-(f4+f16+f18))-1.0/3.0*rho_tmp*uy_tmp;
		Nx_z = 0.5*(f5+f16+f15-(f6+f17+f18))-1.0/3.0*rho_tmp*uz_tmp;
		f[i][j][k][1]  = f2 +1.0*3.0*rho_tmp* ux_tmp;
		f[i][j][k][9]  = f8 +1.0/6.0*rho_tmp*(ux_tmp-uy_tmp) + Nx_y;
		f[i][j][k][7]  = f10+1.0/6.0*rho_tmp*(ux_tmp+uy_tmp) - Nx_y;
		f[i][j][k][11] = f14+1.0/6.0*rho_tmp*(ux_tmp+uz_tmp) - Nx_z;
		f[i][j][k][13] = f12+1.0/6.0*rho_tmp*(ux_tmp-uz_tmp) + Nx_z;
	      }

	      //X-Direction_Outlet
	      if(i==IMAX-1){
		rho_tmp = rho_o;
		uy_tmp  = 0.0; //not neccessarily zero
		uz_tmp  = 0.0; //not neccessarily zero
		ux_tmp   =-1.0+1.0/rho_tmp*(f0+f3+f4+f5+f6+f15+f17+f16+f18+2.0*(f1+f7+f9+f11+f13));
		Nx_y = 0.5*(f3+f15+f17-(f4+f16+f18))-1.0/3.0*rho_tmp*uy_tmp;
		Nx_z = 0.5*(f5+f16+f15-(f6+f17+f18))-1.0/3.0*rho_tmp*uz_tmp;
		f[i][j][k][2]  = f1 -1.0*3.0*rho_tmp* ux_tmp;
		f[i][j][k][8]  = f9 +1.0/6.0*rho_tmp*(-ux_tmp+uy_tmp) - Nx_y;
		f[i][j][k][10] = f7 +1.0/6.0*rho_tmp*(-ux_tmp-uy_tmp) + Nx_y;
		f[i][j][k][14] = f11+1.0/6.0*rho_tmp*(-ux_tmp-uz_tmp) + Nx_z;
		f[i][j][k][12] = f13+1.0/6.0*rho_tmp*(-ux_tmp+uz_tmp) - Nx_z;
		}
	    }
	    if(strcmp(flow_direction,"y")==0){
		//Y-Direction_inlet
	      if(j==0){
		rho_tmp = rho_i;
		ux_tmp  = 0.0; //not neccessarily zero
		uz_tmp  = 0.0; //not neccessarily zero
		uy_tmp  = 1.0-1.0/rho_tmp*(f0+f1+f2+f5+f6+f11+f13+f12+f14+2.0*(f4+f9+f10+f16+f18));
		Ny_x = 0.5*(f1+f11 +f13-(f2+f12+f14))-1.0/3.0*rho_tmp*ux_tmp;
		Ny_z = 0.5*(f5+f11 +f12-(f6+f13+f14))-1.0/3.0*rho_tmp*uz_tmp;
		f[i][j][k][3]  = f4 +1.0*3.0*rho_tmp* uy_tmp;
		f[i][j][k][7]  = f10+1.0/6.0*rho_tmp*(uy_tmp+ux_tmp) - Ny_x;
		f[i][j][k][8]  = f9 +1.0/6.0*rho_tmp*(uy_tmp-ux_tmp) + Ny_x;
		f[i][j][k][15] = f18+1.0/6.0*rho_tmp*(uy_tmp+uz_tmp) - Ny_z;
		f[i][j][k][17] = f16+1.0/6.0*rho_tmp*(uy_tmp-uz_tmp) + Ny_z;
	      }

	      //Y-Direction_Outlet
	      if(j==JMAX-1){
		rho_tmp = rho_o;
		ux_tmp  = 0.0; //not neccessarily zero
		uz_tmp  = 0.0; //not neccessarily zero
		uy_tmp  =-1.0+1.0/rho_tmp*(f0+f1+f2+f5+f6+f9+f10+f13+f14+2.0*(f3+f7+f11+f15+f16));
		Ny_x = 0.5*(f1+f11 +f13-(f2+f12+f14))-1.0/3.0*rho_tmp*ux_tmp;
		Ny_z = 0.5*(f5+f11 +f12-(f6+f13+f14))-1.0/3.0*rho_tmp*uz_tmp;
		f[i][j][k][4]  = f3 -1.0*3.0*rho_tmp*  uy_tmp;
		f[i][j][k][10] = f7 +1.0/6.0*rho_tmp*(-uy_tmp-ux_tmp) + Ny_x;
		f[i][j][k][9]  = f8 +1.0/6.0*rho_tmp*(-uy_tmp+ux_tmp) - Ny_x;
		f[i][j][k][18] = f15+1.0/6.0*rho_tmp*(-uy_tmp-uz_tmp) + Ny_z;
		f[i][j][k][16] = f17+1.0/6.0*rho_tmp*(-uy_tmp+uz_tmp) - Ny_z;
	      }
	    }

	    if(strcmp(flow_direction,"z")==0){
	      //Z-Direction_inlet
	      if(k==0){
		rho_tmp = rho_i;
		ux_tmp  = 0.0; //not neccessarily zero
		uy_tmp  = 0.0; //not neccessarily zero
		uz_tmp  = 1.0-1.0/rho_tmp*(f0+f1+f2+f3+f4+f7+f9+f8+f10+2.0*(f6+f13+f14+f17+f18));
		Nz_x = 0.5*(f1+f7 +f9 -(f2+f8+f10))-1.0/3.0*rho_tmp*ux_tmp;
		Nz_y = 0.5*(f3+f7 +f8 -(f4+f9+f10))-1.0/3.0*rho_tmp*uy_tmp;
		f[i][j][k][5]  = f6 +1.0*3.0*rho_tmp* uz_tmp;
		f[i][j][k][11] = f14+1.0/6.0*rho_tmp*(uz_tmp+ux_tmp) - Nz_x;
		f[i][j][k][12] = f13+1.0/6.0*rho_tmp*(uz_tmp-ux_tmp) + Nz_x;
		f[i][j][k][15] = f18+1.0/6.0*rho_tmp*(uz_tmp+uy_tmp) - Nz_y;
		f[i][j][k][16] = f17+1.0/6.0*rho_tmp*(uz_tmp-uy_tmp) + Nz_y;
	      }

	      //z-Direction_Outlet
	      if(k==KMAX-1){
		rho_tmp = rho_o;
		ux_tmp  = 0.0; //not neccessarily zero
		uy_tmp  = 0.0; //not neccessarily zero
		uz_tmp  =-1.0+1.0/rho_tmp*(f0+f1+f2+f3+f4+f7+f8+f10+f9+2.0*(f5+f11+f12+f15+f16));
		Nz_x = 0.5*(f1+f7 +f9 -(f2+f8+f10))-1.0/3.0*rho_tmp*ux_tmp;
		Nz_y = 0.5*(f3+f7 +f8 -(f4+f9+f10))-1.0/3.0*rho_tmp*uy_tmp;
		f[i][j][k][6]  = f5 -1.0*3.0*rho_tmp* uz_tmp;
		f[i][j][k][13] = f12+1.0/6.0*rho_tmp*(-uz_tmp+ux_tmp) - Nz_x;
		f[i][j][k][14] = f11+1.0/6.0*rho_tmp*(-uz_tmp-ux_tmp) + Nz_x;
		f[i][j][k][17] = f16+1.0/6.0*rho_tmp*(-uz_tmp+uy_tmp) - Nz_y;
		f[i][j][k][18] = f15+1.0/6.0*rho_tmp*(-uz_tmp-uy_tmp) + Nz_y;
	      }
	    }


	  }

	}
}


void macro(void){

  double s_f,s_ux,s_uy,s_uz;

  for(int i=0;i<IMAX;i++)
    for(int j=0;j<JMAX;j++)
      for(int k=0;k<KMAX;k++){

	if(isnode[i][j][k]==0){

	  s_f = 0.0;
	  s_ux = 0.0;
	  s_uy = 0.0;
	  s_uz = 0.0;
	  for(int nv=0;nv<QMAX;nv++){

	    s_f = s_f + f[i][j][k][nv];
	    s_ux = s_ux + f[i][j][k][nv]*eex[nv];
	    s_uy = s_uy + f[i][j][k][nv]*eey[nv];
	    s_uz = s_uz + f[i][j][k][nv]*eez[nv];

	  }

	  rho[i][j][k] = s_f;     //Density
	  ux[i][j][k] = (s_ux+0.5*Fx[i][j][k])/s_f; //Velocity_x
	  uy[i][j][k] = (s_uy+0.5*Fy[i][j][k])/s_f; //Velocity_y
	  uz[i][j][k] = (s_uz+0.5*Fz[i][j][k])/s_f; //Velocity_z

	  press[i][j][k] = cs2*rho[i][j][k]; //Pressure
	}
      }
}

void convergence(void){

  double R2_ux,R2_uy,R2_uz,R2_rho,R2_sum;
  int n_count;
  double R2_uSum;

  R2_ux = 0.0;
  R2_uy = 0.0;
  R2_uz = 0.0;
  R2_sum = 0.0;
  R2_uSum = 0.0;
  R2_rho = 0.0;
  n_count = 0;

  for(int i=0;i<IMAX;i++)
    for(int j=0;j<JMAX;j++)
      for(int k=0;k<KMAX;k++)
	{
	  if(isnode[i][j][k]==0)
	    {
	      R2_ux  += (ux[i][j][k]-ux_old[i][j][k])*(ux[i][j][k]-ux_old[i][j][k]);
	      R2_uy  += (uy[i][j][k]-uy_old[i][j][k])*(uy[i][j][k]-uy_old[i][j][k]);
	      R2_uz  += (uz[i][j][k]-uz_old[i][j][k])*(uz[i][j][k]-uz_old[i][j][k]);
	      R2_rho += (rho[i][j][k]-rho_old[i][j][k])*(rho[i][j][k]-rho_old[i][j][k]);
	      R2_sum += ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k];
	      if(strcmp(flow_direction,"x")==0) R2_uSum += ux[i][j][k]*ux[i][j][k];
	      if(strcmp(flow_direction,"y")==0) R2_uSum += uy[i][j][k]*uy[i][j][k];
	      if(strcmp(flow_direction,"z")==0) R2_uSum += uz[i][j][k]*uz[i][j][k];
	      //	n_count += 1;

	      ux_old[i][j][k] = ux[i][j][k];
	      uy_old[i][j][k] = uy[i][j][k];
	      uz_old[i][j][k] = uz[i][j][k];
	      rho_old[i][j][k] = rho[i][j][k];
	    }

	}

  error_old = error;
  error = sqrt(R2_ux+R2_uy+R2_uz)/sqrt(R2_sum);
  //  error_diff = fabs((error-error_old)/error_old);
  error_diff = fabs((error-error_old));

  printf("errorU: %e |U|: %f\n", error, sqrt(R2_uSum));


}

void results(){

  double allFlow,outFlow,outFlux,outPress,inPress;
  int inCount,outCount,allCount;
  char strK_x[128],strIn[128];



  if(strcmp(vel_output,"on") == 0)
    {
      fpo_vel = fopen("vel_results.dat","w");
      fprintf(fpo_vel,"TITLE= \"3D FLOW\"\n");
      fprintf(fpo_vel,"VARAIBLES = \"I\",  \"J\",  \"K\",  \"VELOCITY_X\"  \"VELOCITY_Y\"  \"VELOCITY_Z\"  \"PRESS\"  \"SOLID\"\n");
      fprintf(fpo_vel,"ZONE T=\"  %d\"  I=   %d,   J=   %d,   K=   %d ",ntime, KMAX-2,JMAX-2,IMAX-10);
      fprintf(fpo_vel,"  F=point \n");
    }


  allFlow = 0.0;
  outFlow = 0.0;
  inCount = 0;
  outCount = 0;
  allCount = 0;
  inPress = 0.0;
  outPress = 0.0;

  if(strcmp(flow_direction,"x")==0){
    for(int i=exLength;i<IMAX-exLength;i++)
      for(int j=1;j<JMAX-1;j++)
	for(int k=1;k<KMAX-1;k++)
	  {
	    if(isnode[i][j][k]==0){
	      allFlow += ux[i][j][k];
	      allCount ++;
	      if(i==exLength){
		inPress += press[i][j][k];
		inCount ++;
	      }
	      if(i==IMAX-exLength-1){
		outPress += press[i][j][k];
		outFlow += ux[i][j][k];
		outFlux += rho[i][j][k]*ux[i][j][k];
		outCount ++;
	      }
	    }
	    else{
	      ux[i][j][k] = 0.0;
	      uy[i][j][k] = 0.0;
	      uz[i][j][k] = 0.0;
	      rho[i][j][k] = 0.0;
	      press[i][j][k] = 0.0;
	    }

	    if(strcmp(vel_output,"on") == 0)
	      {
		fprintf(fpo_vel,"%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%f\t%d\n",i,j,k,ux[i][j][k],uy[i][j][k],uz[i][j][k],press[i][j][k],isnode[i][j][k]);
	      }
	  }
  }

  if(strcmp(flow_direction,"y")==0){
    for(int i=1;i<IMAX-1;i++)
      for(int j=exLength;j<JMAX-exLength;j++)
	for(int k=1;k<KMAX-1;k++)
	  {
	    if(isnode[i][j][k]==0){
	      allFlow += uy[i][j][k];
	      allCount ++;
	      if(j==exLength){
		inPress += press[i][j][k];
		inCount ++;
	      }
	      if(j==JMAX-exLength-1){
		outPress += press[i][j][k];
		outFlow += uy[i][j][k];
		outFlux += rho[i][j][k]*uy[i][j][k];
		outCount ++;
	      }
	    }
	    else{
	      ux[i][j][k] = 0.0;
	      uy[i][j][k] = 0.0;
	      uz[i][j][k] = 0.0;
	      rho[i][j][k] = 0.0;
	      press[i][j][k] = 0.0;
	    }

	    if(strcmp(vel_output,"on") == 0)
	      {
		fprintf(fpo_vel,"%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%f\t%d\n",i,j,k,ux[i][j][k],uy[i][j][k],uz[i][j][k],press[i][j][k],isnode[i][j][k]);
	      }
	  }
  }

  if(strcmp(flow_direction,"z")==0){
    for(int i=1;i<IMAX-1;i++)
      for(int j=1;j<JMAX-1;j++)
	for(int k=exLength;k<KMAX-exLength;k++)
	  {
	    if(isnode[i][j][k]==0){
	      allFlow += uz[i][j][k];
	      allCount ++;
	      if(k==exLength){
		inPress += press[i][j][k];
		inCount ++;
	      }
	      if(k==KMAX-exLength-1){
		outPress += press[i][j][k];
		outFlow += uz[i][j][k];
		outFlux += rho[i][j][k]*uz[i][j][k];
		outCount ++;
	      }
	    }
	    else{
	      ux[i][j][k] = 0.0;
	      uy[i][j][k] = 0.0;
	      uz[i][j][k] = 0.0;
	      rho[i][j][k] = 0.0;
	      press[i][j][k] = 0.0;
	    }

	    if(strcmp(vel_output,"on") == 0)
	      {
		fprintf(fpo_vel,"%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%f\t%d\n",i,j,k,ux[i][j][k],uy[i][j][k],uz[i][j][k],press[i][j][k],isnode[i][j][k]);
	      }
	  }
  }



  outFlow = outFlow/outCount;
  allFlow = allFlow/allCount;
  inPress = inPress/inCount;
  outPress = outPress/outCount;


  double gradP,deltaP,crossArea;

    if(strcmp(PBC_switch,"on")==0){
      if(strcmp(flow_direction,"x")==0){
	gradP = (inPress-outPress)/(double)(IMAX-exLength2-1);
	crossArea = (double)(JMAX-2)*(KMAX-2);
      }
      if(strcmp(flow_direction,"y")==0){
	gradP = (inPress-outPress)/(double)(JMAX-exLength2-1);
	crossArea = (double)(IMAX-2)*(KMAX-2);
      }
      if(strcmp(flow_direction,"z")==0){
	gradP = (inPress-outPress)/(double)(KMAX-exLength2-1);
	crossArea = (double)(IMAX-2)*(JMAX-2);
      }
    }else{
      if(strcmp(flow_direction,"x")==0){
	gradP = fextx;
	deltaP = gradP*(double)(IMAX-1);
	crossArea = (double)(JMAX-2)*(KMAX-2);
      }
      if(strcmp(flow_direction,"y")==0){
	gradP = fexty;
	deltaP = gradP*(double)(JMAX-1);
	crossArea = (double)(IMAX-2)*(KMAX-2);
      }
      if(strcmp(flow_direction,"z")==0){
	gradP = fextz;
	deltaP = gradP*(double)(KMAX-1);
	crossArea = (double)(IMAX-2)*(JMAX-2);
      }
    }



  double K_lbm = outFlow*(double)(outCount/crossArea)*nu/gradP;

  printf("allFlow is %.8f\n",allFlow);
  printf("outFlow is %.8f\n",outFlow);

  printf("PoreCount@all is %d\n",allCount);
  printf("PoreCount@out is %d\n",outCount);
  printf("crossArea@out is %f\n",crossArea);

  ///////////////Result Summary/////////////////

  if(strcmp(vel_output,"on") == 0){

    fpo_sum = fopen("Summary.txt","w");
    fprintf(fpo_sum,"------------------LBM UNITS---------------\n");
    fprintf(fpo_sum,"Lx(lbm) = %i\n",IMAX);
    fprintf(fpo_sum,"Ly(lbm) = %i\n",JMAX);
    fprintf(fpo_sum,"Lz(lbm) = %i\n",KMAX);
    fprintf(fpo_sum,"tTime = %i\n",ntime);
    fprintf(fpo_sum,"Porosity = %f\n",(double)allCount/((IMAX-10)*(JMAX-2)*(KMAX-2)));
    fprintf(fpo_sum,"kinemVis(lbm) = %.2f\n",nu);
    fprintf(fpo_sum,"avgVel(lbm) allPores = %.8f\n",allFlow);
    fprintf(fpo_sum,"avgVel(lbm) outLet = %.8f\n",outFlow);
    fprintf(fpo_sum,"Flux(lbm) outLet = %.8f\n",outFlux);
    fprintf(fpo_sum,"gradP(lbm) outLet = %.8f\n",gradP);
    fprintf(fpo_sum,"deltaP(lbm) outLet = %.8f\n",deltaP);
    fprintf(fpo_sum,"Kx(lbm) = %f\n",K_lbm);
    fprintf(fpo_sum,"Ky(lbm) = %f\n",0.0);
    fprintf(fpo_sum,"Kz(lbm) = %f\n",0.0);

    fprintf(fpo_sum,"\n\n----------------PHYSICAL UNITS---------------\n");
    fprintf(fpo_sum,"pixelSize (m) = %e\n",pixSize);
    fprintf(fpo_sum,"Lx(m) = %e\n",(float)IMAX*pixSize);
    fprintf(fpo_sum,"Ly(m) = %e\n",(float)JMAX*pixSize);
    fprintf(fpo_sum,"Lz(m) = %e\n",(float)KMAX*pixSize);
    fprintf(fpo_sum,"Porosity = %f\n",(double)allCount/((IMAX-10)*(JMAX-2)*(KMAX-2)));
    fprintf(fpo_sum,"Kx(m2) = %e\n", (K_lbm*pixSize*pixSize) );
    fprintf(fpo_sum,"Kx(mDarcy) = %e\n", (K_lbm*pixSize*pixSize)/0.9869233e-12*1.0e3);

    //  fprintf(fpo_sum,"%e\n", (K_lbm*pixSize*pixSize)/0.9869233e-12*1.0e3);
    fclose(fpo_sum);
  }

  fpo_K = fopen(outName,"w");
  fprintf(fpo_K,"{\n");
  fprintf(fpo_K,"\"Porosity\":");
  if(strcmp(flow_direction,"x")==0)
    fprintf(fpo_K,"%f,\n", (double)allCount/(double)((IMAX-exLength2)*(JMAX-2)*(KMAX-2)) );
  if(strcmp(flow_direction,"y")==0)
    fprintf(fpo_K,"%f,\n", (double)allCount/(double)((IMAX-2)*(JMAX-exLength2)*(KMAX-2)) );
  if(strcmp(flow_direction,"z")==0)
    fprintf(fpo_K,"%f,\n", (double)allCount/(double)((IMAX-2)*(JMAX-2)*(KMAX-exLength2)) );
  //  fprintf(fpo_K,"%f,\n", (double)allCount/((IMAX-10)*(JMAX-2)*(KMAX-2)) );
  fprintf(fpo_K,"\"K_lbm\":");
  fprintf(fpo_K,"%e\n", K_lbm);
  fprintf(fpo_K,"}");
  fclose(fpo_K);

}



void free_memory(void){

  int i,j,k;

  for(i=0;i<IMAX;i++)
    for(j=0;j<JMAX;j++)
      for(k=0;k<KMAX;k++)
	{
	  free(f[i][j][k]);
	  free(f_coll[i][j][k]);
	  f[i][j][k] = NULL;
	  f_coll[i][j][k] = NULL;
	}

  for(i=0;i<IMAX;i++){
    for(j=0;j<JMAX;j++)
      {
	free(f[i][j]);
	free(f_coll[i][j]);
	free(ux[i][j]);
	free(uy[i][j]);
	free(uz[i][j]);
	free(rho[i][j]);
	free(press[i][j]);
	free(ux_old[i][j]);
	free(uy_old[i][j]);
	free(uz_old[i][j]);
	free(rho_old[i][j]);
	free(Fx[i][j]);
	free(Fy[i][j]);
	free(Fz[i][j]);
	free(isnode[i][j]);

	f[i][j] = NULL;
	f_coll[i][j] = NULL;
	ux[i][j] = NULL;
	uy[i][j] = NULL;
	uz[i][j] = NULL;
	rho[i][j] = NULL;
	press[i][j] = NULL;
	ux_old[i][j] = NULL;
	uy_old[i][j] = NULL;
	uz_old[i][j] = NULL;
	rho_old[i][j] = NULL;
	Fx[i][j] = NULL;
	Fy[i][j] = NULL;
	Fz[i][j] = NULL;
	isnode[i][j] = NULL;
      }
  }

  for(i=0;i<IMAX;i++){
    free(f[i]);
    free(f_coll[i]);
    free(ux[i]);
    free(uy[i]);
    free(uz[i]);
    free(rho[i]);
    free(press[i]);
    free(ux_old[i]);
    free(uy_old[i]);
    free(uz_old[i]);
    free(rho_old[i]);
    free(Fx[i]);
    free(Fy[i]);
    free(Fz[i]);
    free(isnode[i]);

    f[i] = NULL;
    f_coll[i] = NULL;
    ux[i] = NULL;
    uy[i] = NULL;
    uz[i] = NULL;
    rho[i] = NULL;
    press[i] = NULL;
    ux_old[i] = NULL;
    uy_old[i] = NULL;
    uz_old[i] = NULL;
    rho_old[i] = NULL;
    Fx[i] = NULL;
    Fy[i] = NULL;
    Fz[i] = NULL;
    isnode[i] = NULL;

  }

  free(f);
  free(f_coll);
  free(ux);
  free(uy);
  free(uz);
  free(rho);
  free(press);
  free(ux_old);
  free(uy_old);
  free(uz_old);
  free(rho_old);
  free(Fx);
  free(Fy);
  free(Fz);
  free(isnode);

  f = NULL;
  f_coll = NULL;
  ux = NULL;
  uy = NULL;
  uz = NULL;
  rho = NULL;
  press = NULL;
  ux_old = NULL;
  uy_old = NULL;
  uz_old = NULL;
  rho_old = NULL;
  Fx = NULL;
  Fy = NULL;
  Fz = NULL;
  isnode = NULL;
}

void cvtInt(char *str, char *inP){

  sprintf(str, "%s", inP);
}


void collision_MRT(void)
{

  double u_x,u_y,u_z,uu,arho;
  double Jx,Jy,Jz,JdotJ;
  double meq[QMAX]; // D3Q19 MRT
  double F_x,F_y,F_z,FxUx,FyUy,FzUz,FxUy,FxUz,FyUx,FyUz,FzUx,FzUy,FdotU;
  //  double Fmm1,Fmm2;//(SIMPLIFY)
  //  double sv1_6p,sv7_18p,sv7_14m,sv7_18m,sv11_18m,sv3_6p,sv7_18mm,sv7_14mm;
  //  double fv1_6p,fv7_18p,fv7_14m,fv7_18m,fv11_18m,fv3_6p,fv7_18mm,fv7_14mm;

  double el_cm1,ei_cm1,tw_cm9,tw_cm1012m,tw_cm1012p,tw_cm910p,cm9101112p,cm9101112m,cm0el1m,cm0ei1cm2p;
  double cm012m,cm012m1,cm34p,cm56p;

  for(int i=0;i<IMAX;i++)
    for(int j=0;j<JMAX;j++)
      for(int k=0;k<KMAX;k++)
	{

	  if(isnode[i][j][k]==0){

	    Fx[i][j][k] = fextx;
	    Fy[i][j][k] = fexty;
	    Fz[i][j][k] = fextz;

	    F_x = Fx[i][j][k];
	    F_y = Fy[i][j][k];
	    F_z = Fz[i][j][k];

	    u_x = ux[i][j][k];
	    u_y = uy[i][j][k];
	    u_z = uz[i][j][k];

	    uu = u_x*u_x+u_y*u_y+u_z*u_z;

	    FxUx =F_x*u_x;
	    FyUy =F_y*u_y;
	    FzUz =F_z*u_z;

	    FxUy =F_x*u_y;
	    FxUz =F_x*u_z;

	    FyUx =F_y*u_x;
	    FyUz =F_y*u_z;

	    FzUx =F_z*u_x;
	    FzUy =F_z*u_y;

	    FdotU = FxUx+FyUy+FzUz;

	    Sm[0] = 0.0;
	    Sm[1] = 38.0*FdotU;
	    Sm[2] = -11.0*FdotU;
	    Sm[3] = F_x;
	    Sm[4] = -2.0/3.0*F_x;
	    Sm[5] = F_y;
	    Sm[6] = -2.0/3.0*F_y;
	    Sm[7] = F_z;
	    Sm[8] = -2.0/3.0*F_z;
	    //	    Sm[9] = 2.0*(3.0*FxUx-FdotU);
	    Sm[9] = 2.0*(2.0*FxUx-FyUy-FzUz);
	    //	    Sm[10]= -(3.0*FxUx-FdotU);
	    Sm[10]= -(2.0*FxUx-FyUy-FzUz);
	    Sm[11]= 2.0*(FyUy-FzUz);
	    Sm[12]= -(FyUy-FzUz);
	    Sm[13]= (FxUy+FyUx);
	    Sm[14]= (FyUz+FzUy);
	    Sm[15]= (FxUz+FzUx);
	    Sm[16]= 0.0;
	    Sm[17]= 0.0;
	    Sm[18]= 0.0;


	    for(int nv=0;nv<QMAX;nv++){
	      fv[nv] =f[i][j][k][nv];
	    }


	    // Transforming distribution function into moment space
	    m[0] = fv[0]+fv[1]+fv[2]+fv[3]+fv[4]+fv[5]+fv[6]+fv[7]+fv[8]+fv[9]+fv[10]+fv[11]+fv[12]+fv[13]+fv[14]+fv[15]+fv[16]+fv[17]+fv[18];
	    m[1] = -30.0*fv[0]-11.0*(fv[1]+fv[2]+fv[3]+fv[4]+fv[5]+fv[6])+8.0*(fv[7]+fv[8]+fv[9]+fv[10]+fv[11]+fv[12]+fv[13]+fv[14]+fv[15]+fv[16]+fv[17]+fv[18]);
	    m[2] = 12.0*fv[0]-4.0*(fv[1]+fv[2]+fv[3]+fv[4]+fv[5]+fv[6])+(fv[7]+fv[8]+fv[9]+fv[10]+fv[11]+fv[12]+fv[13]+fv[14]+fv[15]+fv[16]+fv[17]+fv[18]);
	    m[3] = (fv[1]-fv[2]) +fv[7]-fv[8]+fv[9]-fv[10]+fv[11]-fv[12]+fv[13]-fv[14];
	    m[4] = -4.0*(fv[1]-fv[2])+fv[7]-fv[8]+fv[9]-fv[10]+fv[11]-fv[12]+fv[13]-fv[14];
	    m[5] = (fv[3]-fv[4])+(fv[7]+fv[8]-fv[9]-fv[10]) +(fv[15]-fv[16]+fv[17]-fv[18]);
	    m[6] = -4.0*(fv[3]-fv[4])+(fv[7]+fv[8]-fv[9]-fv[10])+(fv[15]-fv[16]+fv[17]-fv[18]);
	    m[7] = (fv[5]-fv[6])+(fv[11]+fv[12]-fv[13]-fv[14]+fv[15]+fv[16]-fv[17]-fv[18]);
	    m[8] = -4.0*(fv[5]-fv[6])+(fv[11]+fv[12]-fv[13]-fv[14]+fv[15]+fv[16]-fv[17]-fv[18]);
	    m[9] = 2.0*(fv[1]+fv[2])-(fv[3]+fv[4]+fv[5]+fv[6])+(fv[7]+fv[8]+fv[9]+fv[10]+fv[11]+fv[12]+fv[13]+fv[14])-2.0*(fv[15]+fv[16]+fv[17]+fv[18]);
	    m[10]= -4.0*(fv[1]+fv[2])+2.0*(fv[3]+fv[4]+fv[5]+fv[6])+(fv[7]+fv[8]+fv[9]+fv[10]+fv[11]+fv[12]+fv[13]+fv[14])-2.0*(fv[15]+fv[16]+fv[17]+fv[18]);
	    m[11]= (fv[3]+fv[4]-fv[5]-fv[6])+fv[7]+fv[8]+fv[9]+fv[10]-fv[11]-fv[12]-fv[13]-fv[14];
	    m[12]= -2.0*(fv[3]+fv[4]-fv[5]-fv[6])+fv[7]+fv[8]+fv[9]+fv[10]-fv[11]-fv[12]-fv[13]-fv[14];
	    m[13]= fv[7]-fv[8]-fv[9]+fv[10];
	    m[14]= fv[15]-fv[16]-fv[17]+fv[18];
	    m[15]= fv[11]-fv[12]-fv[13]+fv[14];
	    m[16]= fv[7]-fv[8]+fv[9]-fv[10]-fv[11]+fv[12]-fv[13]+fv[14];
	    m[17]=-fv[7]-fv[8]+fv[9]+fv[10]+fv[15]-fv[16]+fv[17]-fv[18];
	    m[18]= fv[11]+fv[12]-fv[13]-fv[14]-fv[15]-fv[16]+fv[17]+fv[18];


	    arho = 1.0/rho[i][j][k];
	    Jx  = rho[i][j][k]*u_x;
	    Jy  = rho[i][j][k]*u_y;
	    Jz  = rho[i][j][k]*u_z;
	    JdotJ= (Jx*Jx+Jy*Jy+Jz*Jz);

	    /*
 	    // Equilibrium distribution in moment space
	    meq[0] = rho[i][j][k];                         //rho^eq //density
	    meq[1] = -11.0*rho[i][j][k]+19.0*JdotJ*arho;   //e^eq   //kinetic energy
	    meq[2] = 3.0*rho[i][j][k]-11.0/2.0*JdotJ*arho; //e2^eq  //energy square
	    meq[3] = Jx;                                   //Jx^eq  //x-momentum flux
	    meq[4] = -2.0/3.0*Jx;                          //qx^eq  //x-energy flux
	    meq[5] = Jy;                                   //Jy^eq  //y-momentum flux
	    meq[6] = -2.0/3.0*Jy;                          //qy^eq  //y-energy flux
	    meq[7] = Jz;                                   //Jz^eq  //z-momentum flux
	    meq[8] = -2.0/3.0*Jz;                          //qz^eq  //z-energy flux
	    meq[9] = (3.0*Jx*Jx-JdotJ)*arho;               //pxx^eq //normal viscous stress
	    meq[10]= -0.5*meq[9];                          //txx^eq //-1/2*pxx
	    meq[11]= (Jy*Jy-Jz*Jz)*arho;                   //pww^eq //pyy-pzz
	    meq[12]= -0.5*meq[11];                         //tww^eq //-1/2*txx
	    meq[13]= Jx*Jy*arho;                           //pxy^eq //shear viscous stress
	    meq[14]= Jy*Jz*arho;                           //pyz^eq //shear viscous stress
	    meq[15]= Jx*Jz*arho;                           //pxz^eq //shear viscous stress
	    meq[16]= 0.0;                                  //mx^eq  //3rd-rank tensor
	    meq[17]= 0.0;                                  //my^eq  //3rd-rank tensor
	    meq[18]= 0.0;                                  //mz^eq  //3rd-rank tensor
	    */

	    for(int nv=0;nv<QMAX;nv++){

	      // 19x19 Matrix transforming distribution function into moment space
	      // m[nv] = FtoM(fv[0], fv[1], fv[2], fv[3], fv[4], fv[5], fv[6],fv[7], fv[8], fv[9], fv[10],fv[11],fv[12],fv[13],fv[14],fv[15],fv[16],fv[17],fv[18],nv);

	      // Equilibrium distribution in moment space
	      meq[nv] = Meq(rho[i][j][k],Jx,Jy,Jz,JdotJ,nv);

	      // Relaxation with forcing on in moment space
	      cm[nv]  =-Gamhat[nv]* (m[nv]-meq[nv])+(SCAL[nv]-0.5*Gamhat[nv])*Sm[nv];

	    }


	    double cm0,cm1,cm2,cm3,cm4,cm5,cm6,cm7,cm8,cm9,cm10,cm11,cm12,cm13,cm14,cm15,cm16,cm17,cm18;
	    cm0 = cm[0];
	    cm1 = cm[1];
	    cm2 = cm[2];
	    cm3 = cm[3];
	    cm4 = cm[4];
	    cm5 = cm[5];
	    cm6 = cm[6];
	    cm7 = cm[7];
	    cm8 = cm[8];
	    cm9 = cm[9];
	    cm10 = cm[10];
	    cm11 = cm[11];
	    cm12 = cm[12];
	    cm13 = cm[13];
	    cm14 = cm[14];
	    cm15 = cm[15];
	    cm16 = cm[16];
	    cm17 = cm[17];
	    cm18 = cm[18];

	    cv[0]= cm0-30.0*cm1+12.0*cm2;
	    cv[1]= cm0-11.0*cm1-4.0*(cm2+cm4+cm10)+cm3+2.0*cm9;
	    cv[2]= cm0-11.0*cm1-4.0*(cm2-cm4+cm10)-cm3+2.0*cm9;
	    cv[3]= cm0-11.0*cm1-4.0*(cm2+cm6)+( cm5-cm9+cm11)+2.0*(cm10-cm12);
	    cv[4]= cm0-11.0*cm1-4.0*(cm2-cm6)+(-cm5-cm9+cm11)+2.0*(cm10-cm12);
	    cv[5]= cm0-11.0*cm1-4.0*(cm2+cm8)+( cm7-cm9-cm11)+2.0*(cm10+cm12);
	    cv[6]= cm0-11.0*cm1-4.0*(cm2-cm8)+(-cm7-cm9-cm11)+2.0*(cm10+cm12);
	    cv[7]= cm0+8.0*cm1+cm2+cm3+cm4+cm5+cm6+cm9+cm10+cm11+cm12+cm13+cm16-cm17;
	    cv[8]= cm0+8.0*cm1+cm2-cm3-cm4+cm5+cm6+cm9+cm10+cm11+cm12-cm13-cm16-cm17;
	    cv[9]= cm0+8.0*cm1+cm2+cm3+cm4-cm5-cm6+cm9+cm10+cm11+cm12-cm13+cm16+cm17;
	    cv[10]=cm0+8.0*cm1+cm2-cm3-cm4-cm5-cm6+cm9+cm10+cm11+cm12+cm13-cm16+cm17;
	    cv[11]=cm0+8.0*cm1+cm2+cm3+cm4+cm7+cm8+cm9+cm10-cm11-cm12+cm15-cm16+cm18;
	    cv[12]=cm0+8.0*cm1+cm2-cm3-cm4+cm7+cm8+cm9+cm10-cm11-cm12-cm15+cm16+cm18;
	    cv[13]=cm0+8.0*cm1+cm2+cm3+cm4-cm7-cm8+cm9+cm10-cm11-cm12-cm15-cm16-cm18;
	    cv[14]=cm0+8.0*cm1+cm2-cm3-cm4-cm7-cm8+cm9+cm10-cm11-cm12+cm15+cm16-cm18;
	    cv[15]=cm0+8.0*cm1+cm2+cm5+cm6+cm7+cm8-2.0*(cm9+cm10)+cm14+cm17-cm18;
	    cv[16]=cm0+8.0*cm1+cm2-cm5-cm6+cm7+cm8-2.0*(cm9+cm10)-cm14-cm17-cm18;
	    cv[17]=cm0+8.0*cm1+cm2+cm5+cm6-cm7-cm8-2.0*(cm9+cm10)-cm14+cm17+cm18;
	    cv[18]=cm0+8.0*cm1+cm2-cm5-cm6-cm7-cm8-2.0*(cm9+cm10)+cm14-cm17+cm18;


	    // Post-collision distribution function in velocity space
	    for(int nv=0;nv<QMAX;nv++)
	      {
		//Tranforming the relaxed moments back to velocity space
		//		cv[nv] = MtoF(cm[0], cm[1], cm[2], cm[3], cm[4], cm[5], cm[6],cm[7], cm[8], cm[9], cm[10],cm[11],cm[12],cm[13],cm[14],cm[15],cm[16],cm[17],cm[18],nv);

		//Adding the changes of distribution function due to collision and forcing
		f_coll[i][j][k][nv]=fv[nv]+cv[nv]; // D3Q19 MRT
	      }

	  }
	}

}

double FtoM(double f0,double f1,double f2,double f3,double f4,double f5,double f6,double f7,double f8,double f9,double f10,double f11,double f12,double f13,double f14,double f15,double f16,double f17,double f18,int nq)
{
  double moment;

  switch(nq)
    {

    case 0: moment = f0+f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18;
      break;
    case 1: moment = -30.0*f0-11.0*(f1+f2+f3+f4+f5+f6)+8.0*(f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18);
      break;
    case 2: moment = 12.0*f0-4.0*(f1+f2+f3+f4+f5+f6)+(f7+f8+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18);
      break;
    case 3: moment = (f1-f2) +f7-f8+f9-f10+f11-f12+f13-f14;
      break;
    case 4: moment = -4.0*(f1-f2)+f7-f8+f9-f10+f11-f12+f13-f14;
      break;
    case 5: moment = (f3-f4)+(f7+f8-f9-f10) +(f15-f16+f17-f18);
      break;
    case 6: moment = -4.0*(f3-f4)+(f7+f8-f9-f10)+(f15-f16+f17-f18);
      break;
    case 7: moment = (f5-f6)+(f11+f12-f13-f14+f15+f16-f17-f18);
      break;
    case 8: moment = -4.0*(f5-f6)+(f11+f12-f13-f14+f15+f16-f17-f18);
      break;
    case 9: moment = 2.0*(f1+f2)-(f3+f4+f5+f6)+(f7+f8+f9+f10+f11+f12+f13+f14)-2.0*(f15+f16+f17+f18);
      break;
    case 10:moment = -4.0*(f1+f2)+2.0*(f3+f4+f5+f6)+(f7+f8+f9+f10+f11+f12+f13+f14)-2.0*(f15+f16+f17+f18);
      break;
    case 11:moment = (f3+f4-f5-f6)+f7+f8+f9+f10-f11-f12-f13-f14;
      break;
    case 12:moment = -2.0*(f3+f4-f5-f6)+f7+f8+f9+f10-f11-f12-f13-f14;
      break;
    case 13:moment = f7-f8-f9+f10;
      break;
    case 14:moment = f15-f16-f17+f18;
      break;
    case 15:moment = f11-f12-f13+f14;
      break;
    case 16:moment = f7-f8+f9-f10-f11+f12-f13+f14;
      break;
    case 17:moment =-f7-f8+f9+f10+f15-f16+f17-f18;
      break;
    case 18:moment = f11+f12-f13-f14-f15-f16+f17+f18;
      break;
    default:moment = 0.0;
      break;
    }
  return moment;
}

double Meq(double RHO,double JX,double JY,double JZ,double JDOTJ,int nq)
{

  double MEQ;

  switch(nq)
    {
    case 0: MEQ = RHO;   //rho^eq //density
      break;
    case 1: MEQ = -11.0*RHO+19.0*JDOTJ/RHO;   //e^eq   //kinetic energy
      break;
    case 2: MEQ = 3.0*RHO-11.0/2.0*JDOTJ/RHO; //e2^eq  //energy square
      break;
    case 3: MEQ = JX;           //Jx^eq  //x-momentum flux
      break;
    case 4: MEQ = -2.0/3.0*JX;  //qx^eq  //x-energy flux
      break;
    case 5: MEQ = JY;           //Jy^eq  //y-momentum flux
      break;
    case 6: MEQ = -2.0/3.0*JY;  //qy^eq  //y-energy flux
      break;
    case 7: MEQ = JZ;           //Jz^eq  //z-momentum flux
      break;
    case 8: MEQ = -2.0/3.0*JZ;  //qz^eq  //z-energy flux
      break;
    case 9: MEQ = (3.0*JX*JX-JDOTJ)/RHO;      //pxx^eq //normal viscous stress
      break;
    case 10:MEQ = -0.5*(3.0*JX*JX-JDOTJ)/RHO; //txx^eq //-1/2*pxx
      break;
    case 11:MEQ = (JY*JY-JZ*JZ)/RHO;       //pww^eq //pyy-pzz
      break;
    case 12:MEQ = -0.5*(JY*JY-JZ*JZ)/RHO;  //tww^eq //-1/2*txx
      break;
    case 13:MEQ = JX*JY/RHO;  //pxy^eq //shear viscous stress
      break;
    case 14:MEQ = JY*JZ/RHO;  //pyz^eq //shear viscous stress
      break;
    case 15:MEQ = JX*JZ/RHO;  //pxz^eq //shear viscous stress
      break;
    case 16:MEQ = 0.0;  //mx^eq  //3rd-rank tensor
      break;
    case 17:MEQ = 0.0;  //my^eq  //3rd-rank tensor
      break;
    case 18:MEQ = 0.0;  //mz^eq  //3rd-rank tensor
      break;
    default:MEQ = 0.0;
      break;
    }

  return MEQ;


}

double MtoF(double cm0,double cm1,double cm2,double cm3,double cm4,double cm5,double cm6,double cm7,double cm8,double cm9,double cm10,double cm11,double cm12,double cm13,double cm14,double cm15,double cm16,double cm17,double cm18,int nq)
{
  double pdf;

  switch(nq)
    {

    case 0: pdf = cm0-30.0*cm1+12.0*cm2;
      break;
    case 1: pdf = cm0-11.0*cm1-4.0*(cm2+cm4+cm10)+cm3+2.0*cm9;
      break;
    case 2: pdf = cm0-11.0*cm1-4.0*(cm2-cm4+cm10)-cm3+2.0*cm9;
      break;
    case 3: pdf = cm0-11.0*cm1-4.0*(cm2+cm6)+( cm5-cm9+cm11)+2.0*(cm10-cm12);
      break;
    case 4: pdf = cm0-11.0*cm1-4.0*(cm2-cm6)+(-cm5-cm9+cm11)+2.0*(cm10-cm12);
      break;
    case 5: pdf = cm0-11.0*cm1-4.0*(cm2+cm8)+( cm7-cm9-cm11)+2.0*(cm10+cm12);
      break;
    case 6: pdf = cm0-11.0*cm1-4.0*(cm2-cm8)+(-cm7-cm9-cm11)+2.0*(cm10+cm12);
       break;
    case 7: pdf = cm0+8.0*cm1+cm2+cm3+cm4+cm5+cm6+cm9+cm10+cm11+cm12+cm13+cm16-cm17;
      break;
    case 8: pdf = cm0+8.0*cm1+cm2-cm3-cm4+cm5+cm6+cm9+cm10+cm11+cm12-cm13-cm16-cm17;
      break;
    case 9: pdf = cm0+8.0*cm1+cm2+cm3+cm4-cm5-cm6+cm9+cm10+cm11+cm12-cm13+cm16+cm17;
      break;
    case 10:pdf = cm0+8.0*cm1+cm2-cm3-cm4-cm5-cm6+cm9+cm10+cm11+cm12+cm13-cm16+cm17;
      break;
    case 11:pdf = cm0+8.0*cm1+cm2+cm3+cm4+cm7+cm8+cm9+cm10-cm11-cm12+cm15-cm16+cm18;
      break;
    case 12:pdf = cm0+8.0*cm1+cm2-cm3-cm4+cm7+cm8+cm9+cm10-cm11-cm12-cm15+cm16+cm18;
      break;
    case 13:pdf = cm0+8.0*cm1+cm2+cm3+cm4-cm7-cm8+cm9+cm10-cm11-cm12-cm15-cm16-cm18;
      break;
    case 14:pdf = cm0+8.0*cm1+cm2-cm3-cm4-cm7-cm8+cm9+cm10-cm11-cm12+cm15+cm16-cm18;
      break;
    case 15:pdf = cm0+8.0*cm1+cm2+cm5+cm6+cm7+cm8-2.0*(cm9+cm10)+cm14+cm17-cm18;
      break;
    case 16:pdf = cm0+8.0*cm1+cm2-cm5-cm6+cm7+cm8-2.0*(cm9+cm10)-cm14-cm17-cm18;
      break;
    case 17:pdf = cm0+8.0*cm1+cm2+cm5+cm6-cm7-cm8-2.0*(cm9+cm10)-cm14+cm17+cm18;
      break;
    case 18:pdf = cm0+8.0*cm1+cm2-cm5-cm6-cm7-cm8-2.0*(cm9+cm10)+cm14-cm17+cm18;
      break;
    default:pdf = 0.0;
      break;
    }

  return pdf;
}


// Read files
void readfile(char* filepath, char* fileContent)
{
  FILE *f;
  char c;
  int index;

  f = fopen(filepath, "rt");
  while((c = fgetc(f)) != EOF){
    fileContent[index] = c;
    index++;
  }
  fileContent[index] = '\0';
}


int parseJSON(char *filepath, void callback(char *, char*)){

  char JSON_STRING[BUFFER_SIZE];

  char value[1024];
  char key[1024];

  //read json file
  readfile(filepath, JSON_STRING);

  int i;
  int r;

  jsmn_parser p;
  jsmntok_t t[MAX_TOKEN_COUNT];

  jsmn_init(&p);

  r = jsmn_parse(&p, JSON_STRING, strlen(JSON_STRING), t, sizeof(t)/(sizeof(t[0])));

  if (r < 0) {
    printf("Failed to parse JSON: %d\n", r);
    return 1;
  }

  /* Assume the top-level element is an object */
  if (r < 1 || t[0].type != JSMN_OBJECT) {
    printf("Object expected\n");
    return 1;
  }

  for (i = 1; i < r; i++){

    jsmntok_t json_value = t[i+1];
    jsmntok_t json_key = t[i];


    int string_length = json_value.end - json_value.start;
    int key_length = json_key.end - json_key.start;

    int idx;

    for (idx = 0; idx < string_length; idx++){
      value[idx] = JSON_STRING[json_value.start + idx ];
      //           printf("%i\n",value[idx]);
    }

    for (idx = 0; idx < key_length; idx++){
      key[idx] = JSON_STRING[json_key.start + idx];
      //           printf("%i\n",key[idx]);
    }

    value[string_length] = '\0';
    key[key_length] = '\0';

    callback(key, value);

    if(strcmp(key,"flow_direction") == 0) strcpy(flow_direction,value);
    if(strcmp(key,"imax") == 0) IMAX = atoi(value);
    if(strcmp(key,"jmax") == 0) JMAX = atoi(value);
    if(strcmp(key,"kmax") == 0) KMAX = atoi(value);
    if(strcmp(key,"tmax") == 0) TMAX = atoi(value);
    if(strcmp(key,"relaxation_time") == 0) tau = atof(value);
    if(strcmp(key,"viscosity") == 0) nu = atof(value);
    if(strcmp(key,"external_force") == 0) force = atof(value);
    if(strcmp(key,"pressure_gradient_x") == 0) fextx = atof(value);
    if(strcmp(key,"pressure_gradient_y") == 0) fexty = atof(value);
    if(strcmp(key,"pressure_gradient_z") == 0) fextz = atof(value);
    if(strcmp(key,"convergence_tolerance") == 0) tolerance = atof(value);
    if(strcmp(key,"is_x_periodic") == 0) strcpy(periodic_x,value);
    if(strcmp(key,"is_y_periodic") == 0) strcpy(periodic_y,value);
    if(strcmp(key,"is_z_periodic") == 0) strcpy(periodic_z,value);
    if(strcmp(key,"output_velocity_fields") == 0) strcpy(vel_output,value);
    if(strcmp(key,"pixel_resolution(meter)") == 0) pixSize = atof(value);
    if(strcmp(key,"collision_model") == 0) strcpy(collision_model,value);
    if(strcmp(key,"input_data_format") == 0) strcpy(dataFormat,value);
    if(strcmp(key,"input_file_name") == 0) strcpy(input_file_name,value);
    if(strcmp(key,"output_file_name") == 0) strcpy(output_file_name,value);
    if(strcmp(key,"pressure_BC") == 0) strcpy(PBC_switch,value);
    if(strcmp(key,"inlet_density") == 0) rho_i=atof(value);
    if(strcmp(key,"outlet_density") == 0) rho_o=atof(value);

    i++;

  }

  return 0;
}

// Only prints the key and value
void mycallback(char *key, char* value){

  //  printf("%s : %s\n", key, value);

}
