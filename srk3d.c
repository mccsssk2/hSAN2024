/*
2 December 2019.

Morph of TH shemer/3Dsrk.c for human SAN 3D.

I think reason for past divergence was that shape revision I made.
*/

static char help[] = "srk3d.c, Human SAN calculation. \n";

// dimensions as came out of hSAN.m
/*
#define usr_MX 	    218
#define usr_MY       223
#define usr_MZ	    231
*/

// revised real geometry (not idealised).
#define usr_MX 	    216
#define usr_MY      223
#define usr_MZ	    220

// standard headers, see if there is anything ttnp specific (probably not).
// space step is 0.175 mm.
// #include "ttnp_2D.h"
#include "../../san.h"
#include"../../f_SAN.c"
#include "../../f_atrial.c"
extern PetscErrorCode FormIFunction(TS,PetscReal,Vec,Vec,Vec,void*);

int main(int argc,char **argv)
{
	  TS             	ts;                   			/* nonlinear solver 			*/
	  Vec            	u, r;                  		/* solution, residual vectors 	*/
	  Mat            	J	, Jmf = NULL  ;   	/* Jacobian matrices 			*/
//	  PetscInt       maxsteps = 10000;      	/* iterations for convergence 	*/
	  DM             	da;					// a short form of user.da in the main driver function.
	  AppCtx3D    user;              			/* user-defined work context 	*/
	  SNES           	snes;
	// sk variables.
	  PetscInt       mybase_x, mybase_y, mybase_z, mysize_x, mysize_y, mysize_z, usr_i, usr_j, usr_k, time_int, total_time_int, ierr;
	  PetscScalar    ***u_localptr;
	  char 		str[1000]; // so you do not keep creating and destroying this.
	  double   	****usr_u0_loc; // state variable array independent of Petsc and Sundials.
	  realtype       usr_t, usr_tout;
          PetscInt       file_Counter = 0;
          int 		intx, inty, intz;
//         FILE       	*geometry, *output, *spiralFile;
         FILE       	*geometry, *output;
         int 			int0, int1, int2, int3, int4, int5, int6, neq;
         double 		stim_time = 0.0;
         double 		kko, qvalue;
// double in1, in2, in3, in4, in5, in6, in7, in8, in9, in10, in11, in12, in13, in14, in15, in16, in17, in18, in19, in20, in21, in22, in23, in24, in25, in26, in27, in28, in29, in30, in31, in32, in33, in34, in35, in36, in37, in38, in39;
double BCL, usr_time;

	  PetscInitialize(&argc,&argv,(char*)0,help);
	  
	time_int = 0; total_time_int = (int)(10.0 / DELTAT); // 10 s run.
	stim_time = 0.0; usr_time = 0.0;
	for(time_int = 0; time_int < total_time_int; time_int++){ 

	     /* Initialize user application context */
	      user.da           = NULL;
DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,usr_MX,usr_MY,usr_MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,NULL,&da);  /* this decides that u is non-ghosted 1st order star stencil vector. */
	  ierr = DMSetUp(da);CHKERRQ(ierr);
	     user.da = da		;
	     DMCreateGlobalVector(da,&u); 
	     VecDuplicate(u,&r)		; // this is the residual vector, which has the same properties as the state vector u. Residual has to be minimised.
	   DMDAGetCorners(da,&mybase_x,&mybase_y,&mybase_z,&mysize_x,&mysize_y,&mysize_z);  // you need this for running the Sundials loops.

BCL = 0.3; // seconds.

	     if(time_int<1){
		     
	usr_u0_loc 		= (double 	****	) calloc(mysize_z, sizeof(double***));
	user.geometry 	= (PetscInt 	****	) calloc(mysize_z, sizeof(PetscInt***));
	user.data 		= (PetscReal 	****	) calloc(mysize_z, sizeof(PetscReal***));		
	user.dist 		= (double 	***	) calloc(mysize_z, sizeof(double**));			
	
	for (usr_k = 0; usr_k < mysize_z; usr_k++){ 
	usr_u0_loc[usr_k] 	= (double 	***	) calloc(mysize_y, sizeof(double**));
	user.geometry[usr_k] 	= (PetscInt 	***	) calloc(mysize_y, sizeof(PetscInt**));				
	user.data[usr_k] 		= (PetscReal 	***	) calloc(mysize_y, sizeof(PetscReal**));			
	user.dist[usr_k] 		= (double**	) calloc(mysize_y, sizeof(double*));				
				/* the usr_mysize_ and usr_mybase_ arrays are not the same on all processors. */
				for (usr_j = 0; usr_j < mysize_y; usr_j++){ 
				usr_u0_loc[usr_k][usr_j] 		= (double	**	) calloc(mysize_x,sizeof(double     *));
				user.geometry[usr_k][usr_j] 	= (PetscInt	**	) calloc(mysize_x,sizeof(PetscInt   *));
				user.data[usr_k][usr_j] 		= (PetscReal	**	) calloc(mysize_x,sizeof(PetscReal*));				
				user.dist[usr_k][usr_j] 		= (double*	) calloc(mysize_x,sizeof(double));								
				
					for (usr_i = 0; usr_i < mysize_x; usr_i++){
						usr_u0_loc[usr_k][usr_j][usr_i] 		= (double       *) calloc(NEQ_SAN,    sizeof(double     ));						
						user.geometry[usr_k][usr_j][usr_i] 	= (PetscInt	*) calloc(NBS3D,    sizeof(PetscInt   ));										
						user.data[usr_k][usr_j][usr_i] 		= (PetscReal	*) calloc(DATA3D, sizeof(PetscReal));										
					}
				}
	}

geometry = fopen("humanSAN2019.geom","r");
while(fscanf(geometry,"%d %d %d %d  %d %d %d  %d %d %d %lf %lf", &intz, &inty, &intx, &int0, &int1, &int2, &int3, &int4, &int5, &int6, &kko, &qvalue)!=EOF){
				usr_k = (PetscInt)intz - mybase_z; usr_j = (PetscInt)inty - mybase_y; usr_i = (PetscInt)intx  - mybase_x;
				if(usr_k>=0&&usr_k<mysize_z&&usr_j>=0&&usr_j<mysize_y&&usr_i>=0&&usr_i<mysize_x){
				user.geometry[usr_k][usr_j][usr_i][0] = (PetscInt)int0; // 7 types of neighbours: me, top, right, bottom, left.
				user.geometry[usr_k][usr_j][usr_i][1] = (PetscInt)int1; 
				user.geometry[usr_k][usr_j][usr_i][2] = (PetscInt)int2; 
				user.geometry[usr_k][usr_j][usr_i][3] = (PetscInt)int3; 
				user.geometry[usr_k][usr_j][usr_i][4] = (PetscInt)int4; 
				user.geometry[usr_k][usr_j][usr_i][5] = (PetscInt)int5; 
				user.geometry[usr_k][usr_j][usr_i][6] = (PetscInt)int6; 

				// model data.
				user.dist[usr_k][usr_j][usr_i] = (PetscReal)kko; // distance.
				if(argc==3){
				user.data[usr_k][usr_j][usr_i][0] = (PetscReal)qvalue*atof(argv[2]);
//				printf("ical is %f\n", atof(argv[2]));
				}
				else{
				user.data[usr_k][usr_j][usr_i][0] = (PetscReal)qvalue;				
//				printf("ical is 1\n");
				}
				

				} // end of if statement
			} // end of while statement

fclose(geometry);

if(1==2){

		sprintf(str,"HumanSANDistanceinSRK3D.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 1.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX, usr_MY, usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS voltage float 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++){
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
		fprintf(output,"%f ",user.dist[usr_k][usr_j][usr_i]);
	fprintf(output,"\n");
	}
		fclose(output);

}		

// printf("read in geomtry.\n");
// exit(0);

	/* Sundials initial conditions                                                             */
		for (usr_k = 0; usr_k < mysize_z; usr_k++)
		for(usr_j  = 0; usr_j  < mysize_y; usr_j++)
		for(usr_i  = 0; usr_i  < mysize_x; usr_i++){

			for(neq  = 0; neq < NEQ_SAN; neq++)  usr_u0_loc[usr_k][usr_j][usr_i][neq] = 0.0;

if(user.geometry[usr_k][usr_j][usr_i][0]==1){ // Either 1 (atrial), 2 (SAN).
   usr_u0_loc[usr_k][usr_j][usr_i][0] /* Ithy_atrial, 0+1) */ =  0.632613;   // Ca_rel (millimolar) (in Ca_handling_by_the_SR)
   usr_u0_loc[usr_k][usr_j][usr_i][1] /* Ithy_atrial, 1+1) */ =  0.649195;   // Ca_up (millimolar) (in Ca_handling_by_the_SR)
   usr_u0_loc[usr_k][usr_j][usr_i][2] /* Ithy_atrial, 2+1) */ =  0.470055;   // F1 (dimensionless) (in Ca_handling_by_the_SR)
   usr_u0_loc[usr_k][usr_j][usr_i][3] /* Ithy_atrial, 3+1) */ =  0.002814;   // F2 (dimensionless) (in Ca_handling_by_the_SR)
   usr_u0_loc[usr_k][usr_j][usr_i][4] /* Ithy_atrial, 4+1) */ =  0.431547;   // O_Calse (dimensionless) (in Ca_handling_by_the_SR)
   usr_u0_loc[usr_k][usr_j][usr_i][5] /* Ithy_atrial, 5+1) */ =  0.001089;   // r (dimensionless) (in Ca_independent_transient_outward_K_current_r_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][6] /* Ithy_atrial, 6+1) */ =  0.948597;   // s (dimensionless) (in Ca_independent_transient_outward_K_current_s_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][7] /* Ithy_atrial, 7+1) */ =  0.000014;   // d_L (dimensionless) (in L_type_Ca_channel_d_L_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][8] /* Ithy_atrial, 8+1) */ =  0.998597;   // f_L1 (dimensionless) (in L_type_Ca_channel_f_L1_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][9] /* Ithy_atrial, 9+1) */ =  0.998586;   // f_L2 (dimensionless) (in L_type_Ca_channel_f_L2_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][10] /* Ithy_atrial, 10+1) */ =  1.815768;   // Ca_c (millimolar) (in cleft_space_ion_concentrations)
   usr_u0_loc[usr_k][usr_j][usr_i][11] /* Ithy_atrial, 11+1) */ =  5.560224;   // K_c (millimolar) (in cleft_space_ion_concentrations)
   usr_u0_loc[usr_k][usr_j][usr_i][12] /* Ithy_atrial, 12+1) */ =  130.022096;   // Na_c (millimolar) (in cleft_space_ion_concentrations)
   usr_u0_loc[usr_k][usr_j][usr_i][13] /* Ithy_atrial, 13+1) */ =  0.004374;   // n (dimensionless) (in delayed_rectifier_K_currents_n_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][14] /* Ithy_atrial, 14+1) */ =  0.000053;   // pa (dimensionless) (in delayed_rectifier_K_currents_pa_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][15] /* Ithy_atrial, 15+1) */ =  1.38222;   // O (dimensionless) (in intracellular_Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][16] /* Ithy_atrial, 16+1) */ =  0.026766;   // O_C (dimensionless) (in intracellular_Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][17] /* Ithy_atrial, 17+1) */ =  0.012922;   // O_TC (dimensionless) (in intracellular_Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][18] /* Ithy_atrial, 18+1) */ =  0.190369;   // O_TMgC (dimensionless) (in intracellular_Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][19] /* Ithy_atrial, 19+1) */ =  0.714463;   // O_TMgMg (dimensionless) (in intracellular_Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][20] /* Ithy_atrial, 20+1) */ =  7.1e-5;   // Ca_d (millimolar) (in intracellular_ion_concentrations)
   usr_u0_loc[usr_k][usr_j][usr_i][21] /* Ithy_atrial, 21+1) */ =  6.5e-5;   // Ca_i (millimolar) (in intracellular_ion_concentrations)
   usr_u0_loc[usr_k][usr_j][usr_i][22] /* Ithy_atrial, 22+1) */ =  129.485991;   // K_i (millimolar) (in intracellular_ion_concentrations)
   usr_u0_loc[usr_k][usr_j][usr_i][23] /* Ithy_atrial, 23+1) */ =  8.516766;   // Na_i (millimolar) (in intracellular_ion_concentrations)
   usr_u0_loc[usr_k][usr_j][usr_i][24] /* Ithy_atrial, 24+1) */ =  -74.031982;   // V (millivolt) (in membrane)
   usr_u0_loc[usr_k][usr_j][usr_i][25] /* Ithy_atrial, 25+1) */ =  0.877202;   // h1 (dimensionless) (in sodium_current_h1_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][26] /* Ithy_atrial, 26+1) */ =  0.873881;   // h2 (dimensionless) (in sodium_current_h2_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][27] /* Ithy_atrial, 27+1) */ =  0.003289;   // m (dimensionless) (in sodium_current_m_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][28] /* Ithy_atrial, 28+1) */ =  0.000367;   // a_ur (dimensionless) (in ultra_rapid_K_current_aur_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][29] /* Ithy_atrial, 29+1) */ =  0.96729;   // i_ur (dimensionless) (in ultra_rapid_K_current_iur_gate)

   usr_u0_loc[usr_k][usr_j][usr_i][30] =  0.0;   // atrial, dummy.
   usr_u0_loc[usr_k][usr_j][usr_i][31] =  0.0;   // atrial, dummy.
   usr_u0_loc[usr_k][usr_j][usr_i][32] =  0.0;   // atrial, dummy.
}
		
if(user.geometry[usr_k][usr_j][usr_i][0]==2){ // Either 1 (atrial), 2 (SAN).
   //---------------------------------------------------------------------------
   // State variables for SAN
   //---------------------------------------------------------------------------
   usr_u0_loc[usr_k][usr_j][usr_i][0] =   4.595622e-10;   // I (dimensionless) (in Ca_SR_release)
   usr_u0_loc[usr_k][usr_j][usr_i][1] =   6.181512e-9;   // O (dimensionless) (in Ca_SR_release)
   usr_u0_loc[usr_k][usr_j][usr_i][2] =   0.9308;   // R1 (dimensionless) (R in Ca_SR_release)
   usr_u0_loc[usr_k][usr_j][usr_i][3] =   0.069199;   // RI (dimensionless) (in Ca_SR_release)
   usr_u0_loc[usr_k][usr_j][usr_i][4] =   0.217311;   // fCMi (dimensionless) (in Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][5] =   0.158521;   // fCMs (dimensionless) (in Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][6] =   0.138975;   // fCQ (dimensionless) (in Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][7] =   0.017929;   // fTC (dimensionless) (in Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][8] =   0.259947;   // fTMC (dimensionless) (in Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][9] =   0.653777;   // fTMM (dimensionless) (in Ca_buffering)
   usr_u0_loc[usr_k][usr_j][usr_i][10] =   0.409551;   // Ca_jsr (millimolar) (in Ca_dynamics)
   usr_u0_loc[usr_k][usr_j][usr_i][11] =   0.435148;   // Ca_nsr (millimolar) (in Ca_dynamics)
   usr_u0_loc[usr_k][usr_j][usr_i][12] =   6.226104e-5;   // Ca_sub (millimolar) (in Ca_dynamics)
   usr_u0_loc[usr_k][usr_j][usr_i][13] =   9.15641e-6;   // Cai (millimolar) (in Ca_dynamics)
   usr_u0_loc[usr_k][usr_j][usr_i][14] =   -47.787168;   // V_ode (millivolt) (in Membrane)
   usr_u0_loc[usr_k][usr_j][usr_i][15] =   5.0;   // Nai_ (millimolar) (in Nai_concentration)
   usr_u0_loc[usr_k][usr_j][usr_i][16] =   0.001921;   // dL (dimensionless) (in i_CaL_dL_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][17] =   0.844449;   // fCa (dimensionless) (in i_CaL_fCa_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][18] =   0.846702;   // fL (dimensionless) (in i_CaL_fL_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][19] =   0.268909;   // dT (dimensionless) (in i_CaT_dT_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][20] =   0.020484;   // fT (dimensionless) (in i_CaT_fT_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][21] =   0.00277;   // a (dimensionless) (in i_KACh_a_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][22] =   0.011068;   // paF (dimensionless) (in i_Kr_pa_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][23] =   0.283185;   // paS (dimensionless) (in i_Kr_pa_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][24] =   0.709051;   // piy (dimensionless) (in i_Kr_pi_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][25] =   0.1162;   // n (dimensionless) (in i_Ks_n_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][26] =   0.011845;   // r_Kur (dimensionless) (in i_Kur_rKur_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][27] =   0.845304;   // s_Kur (dimensionless) (in i_Kur_sKur_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][28] =   0.003058;   // h (dimensionless) (in i_Na_h_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][29] =   0.447724;   // m (dimensionless) (in i_Na_m_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][30] =   0.009508;   // y (dimensionless) (in i_f_y_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][31] =   0.430836;   // q (dimensionless) (in i_to_q_gate)
   usr_u0_loc[usr_k][usr_j][usr_i][32] =   0.014523;   // r (dimensionless) (in i_to_r_gate)
} // end of SAN init.

}  // end of usr_i usr_j usr_k loops.
//		S2 = 0;
	     } // end of time_int < 1 ****************************************************************************************************************************************
// S2 here if you want.

//******************************************************************************************************************************************************************************
	for (usr_k = 0; usr_k < mysize_z; usr_k++)
	for(usr_j = 0; usr_j < mysize_y; usr_j++) // these usr_j and usr_i are for local locations.
	for(usr_i = 0; usr_i < mysize_x; usr_i++){

		if(user.geometry[usr_k][usr_j][usr_i][0]>0){
	/*******************************************************************************/
		    N_Vector y;                    /* The state vector                 */
		    void *cvode_mem;               /* Sundials memory pointer          */
		    UserData data;                 /* Sundials data structure instance */ 
		    /* Create serial vector of length NEQ for I.C. and abstol          */
		    cvode_mem 	= NULL;
		    data      		= NULL;
		    y         		= NULL;
		    
		    if(user.geometry[usr_k][usr_j][usr_i][0]==1) y         	= N_VNew_Serial(NEQA);
		    if(user.geometry[usr_k][usr_j][usr_i][0]==2) y         	= N_VNew_Serial(NEQ_SAN);
		    
		    data      						= (UserData) malloc(sizeof *data);
		    data->Istim 					= 0.0; // pA/pF.
		    data->extracellularPotassium 	= 5.4; // mM.
		    data->qvalue 					= user.data[usr_k][usr_j][usr_i][0];
		data->inamultiplier = atof(argv[1]); // not in the 1.Control simulations.		    

// if(usr_time<=4.0)
// if(usr_time<=-1.0)
// if(usr_time<=4.0)
// pacing is on.
if(usr_time<=4.0&&stim_time>0.0&&stim_time<3*0.006&&(((usr_i+mybase_x-42)*(usr_i+mybase_x-42)+(usr_j+mybase_y-100)*(usr_j+mybase_y-100) +(usr_k+mybase_z-170)*(usr_k+mybase_z-170))<4*4)) data->Istim=-15.0; else data->Istim=0.0;

// pacing is off.
/*
data->Istim = 0.0;
if(usr_time<=-1.0&&stim_time>0.0&&stim_time<3*0.006&&(((usr_i+mybase_x-42)*(usr_i+mybase_x-42)+(usr_j+mybase_y-100)*(usr_j+mybase_y-100) +(usr_k+mybase_z-170)*(usr_k+mybase_z-170))<4*4)) data->Istim = -15.0; else data->Istim = 0.0;
*/

	    /* user specified Sundials data.                */			
	if(user.geometry[usr_k][usr_j][usr_i][0]==1)
		for(neq=0;neq<NEQA;neq++) 			Ith(y,neq+1) = usr_u0_loc[usr_k][usr_j][usr_i][neq];
	if(user.geometry[usr_k][usr_j][usr_i][0]==2)
		for(neq=0;neq<NEQ_SAN;neq++) 			Ith(y,neq+1) = usr_u0_loc[usr_k][usr_j][usr_i][neq];

		    usr_t     = 0; usr_tout = DELTAT;
		    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		    
	    	if(user.geometry[usr_k][usr_j][usr_i][0]==1)
		    CVodeInit(cvode_mem, f_atrial, 0.0, y);
	    	if(user.geometry[usr_k][usr_j][usr_i][0]==2)
		    CVodeInit(cvode_mem, f_SAN, 0.0, y);

	    	if(user.geometry[usr_k][usr_j][usr_i][0]==1)
		    CVodeWFtolerances(cvode_mem, ewt1); // the NEQ has to be two ways here.
	    	if(user.geometry[usr_k][usr_j][usr_i][0]==2)
		    CVodeWFtolerances(cvode_mem, ewt2); // the NEQ has to be two ways here.

	    CVodeSetUserData(cvode_mem, data);
	    CVodeSStolerances(cvode_mem, RTOL, ATOL);
	    	if(user.geometry[usr_k][usr_j][usr_i][0]==1)
			CVDense(cvode_mem, NEQA);
	    	if(user.geometry[usr_k][usr_j][usr_i][0]==2)
			CVDense(cvode_mem, NEQ_SAN);
		    
		    CVodeSetMaxStep(cvode_mem, DELTAT);
		    CVode(cvode_mem, usr_tout, y, &usr_t, CV_NORMAL);
		    
	if(user.geometry[usr_k][usr_j][usr_i][0]==1)		    
		for(neq=0;neq<NEQA;neq++) 													usr_u0_loc[usr_k][usr_j][usr_i][neq] = Ith(y,neq+1);
	if(user.geometry[usr_k][usr_j][usr_i][0]==2)
		for(neq=0;neq<NEQ_SAN;neq++) 												usr_u0_loc[usr_k][usr_j][usr_i][neq] = Ith(y,neq+1);

		    N_VDestroy_Serial(y);
		    CVodeFree(&cvode_mem);
		    free(data);
	/*******************************************************************************/
	} /* end of geometry if.                                                        */
} // end of Sundials loop.

			// Petsc part, solve for 1 time step from t = 0 to t = DELTAT.
			//     Create timestepping solver context
		   	TSCreate(PETSC_COMM_WORLD,&ts)		; 
		   	TSSetProblemType(ts,TS_NONLINEAR)		; 
		   	TSSetType(ts,TSBEULER)					; 
		   	TSSetDM(ts,da)							;
		   	TSSetIFunction(ts,r,FormIFunction,&user)	; 
		   	
//		   	TSSetDuration(ts,MAXSTEPS,DELTAT)		; // this needs to be the time step/output step of my application.   // this is depercated, but it worked on Cedar on April 25, 2021.
			  ierr = TSSetMaxTime(ts,DELTAT);CHKERRQ(ierr); // replaces TSSetDuration, but I dont know where to use MAXSTEPS.
		   	
			TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); // get the solution at the exact time of end. This is new to LHSC SRK*.C codes.

			//     Set initial conditions ********************************************
    	DMDAVecGetArray(da,u,&u_localptr);
	for (usr_k = 0; usr_k < mysize_z; usr_k++) for(usr_j = 0; usr_j < mysize_y; usr_j++) for(usr_i = 0; usr_i < mysize_x; usr_i++){
			 if(user.geometry[usr_k][usr_j][usr_i][0]==1)				
			u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_k][usr_j][usr_i][24];
			 else if(user.geometry[usr_k][usr_j][usr_i][0]==2)
			u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_k][usr_j][usr_i][14];
	}
	
    	DMDAVecRestoreArray(da,u,&u_localptr);
                        //  ******************************************************************************

			   TSSetSolution(ts,u); 
//			   TSSetInitialTimeStep(ts,0.0,DELTAT); // if convergence fails, then you can reduce DELTAT by several fractions.
//			   TSSetInitialTimeStep(ts,0.0,DELTAT); // this is depercated, but it worked on Cedar on April 25, 2021.
			  ierr = TSSetTimeStep(ts,DELTAT);CHKERRQ(ierr); // replaces TSSetInitialTimeStep.

				//   Set Jacobian evaluation routine using colouring.
			   DMSetMatType(da,MATAIJ); 
			   DMCreateMatrix(da,&J); 
			   TSGetSNES(ts,&snes); 
			   MatCreateSNESMF(snes,&Jmf); 
			   SNESSetJacobian(snes,Jmf,J,SNESComputeJacobianDefaultColor,0);  // should I do some ghosting here? Does not seem necessary.
				//   Sets various TS parameters from user options
			   TSSetFromOptions(ts); 
				// get a solution from 0 to DELTAT
			   TSSolve(ts,u); 
				// get the u into my non PetSc array.
			    DMDAVecGetArray(da,u,&u_localptr);
    			    for (usr_k = 0; usr_k < mysize_z; usr_k++) for(usr_j = 0; usr_j < mysize_y; usr_j++) for(usr_i = 0; usr_i < mysize_x; usr_i++){
			 if(user.geometry[usr_k][usr_j][usr_i][0]==1)					                            
    				usr_u0_loc[usr_k][usr_j][usr_i][24] = u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i];
			 else if(user.geometry[usr_k][usr_j][usr_i][0]==2)				                                
    				usr_u0_loc[usr_k][usr_j][usr_i][14] = u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i];                                
                            }
    				DMDAVecRestoreArray(da,u,&u_localptr);

				// do a vts view.
				if(time_int % (int)(0.001/DELTAT) == 0 /* && usr_time >= 3.5 */ ){ // every 1 ms.
/*
					sprintf(str,"my_3d%d.vts", (int)file_Counter);
					PetscViewer viewer;  PetscViewerCreate(PETSC_COMM_WORLD, &viewer); PetscViewerSetType(viewer, PETSCVIEWERVTK);
					PetscViewerFileSetName(viewer, str); VecView(u, viewer); PetscViewerDestroy(&viewer);
*/

		        sprintf(str,"my_3d%d.bin",(int)file_Counter);
			PetscViewer viewer2;
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_WRITE,&viewer2);
			VecView(u,viewer2);
			PetscViewerDestroy(&viewer2);					

			file_Counter++;
				}

			   MatDestroy(&J); 
			   MatDestroy(&Jmf); 
			   VecDestroy(&u); 
			   VecDestroy(&r); 
			   TSDestroy(&ts); 
			   DMDestroy(&da); 

			// track time
			usr_time = DELTAT * time_int; stim_time = stim_time + DELTAT;
			if(stim_time > BCL) stim_time = 0.0;

} // end of time loop.

// necessary.
   PetscFinalize();
  PetscFunctionReturn(0);
} // end of main function.

/* --------------------------------------------------------------------- */
/*
  FormIFunction = Udot - RHSFunction, RHS minus RHSFunction.
*/
PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx)
{
//  PetscErrorCode ierr;
  AppCtx3D         *user=(AppCtx3D*)ctx;
  DM             da   = (DM)user->da;
  PetscInt       i,j, k, xs,ys, zs, xm,ym, zm, gi, gj, gk;

  PetscScalar    u, uxx, uyy, uzz, uxy, ***uarray, ***f, ***udot; // this is 3D.
  Vec            localU; // petsc vectors are always 1D. I need to cast my vectors like that.
// PetscReal ddxx, ddyy, ddzz, ddxy, fx, fy, fz;
 PetscReal      sx, sy, sz;
// PetscReal U0, U1, U2, U3, U4, U5, U6;
// PetscReal nx, ny, nz;
// PetscReal ux, uy, uz;
// bcsx, bcsy, bcsz;

    PetscFunctionBeginUser;
    DMGetLocalVector(da,&localU); 

  /* Scatter ghost points to local vector,using the 2-step process DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be done while messages are in transition. */
   DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);  DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU); 
  /* Get pointers to vector data */
   DMDAVecGetArrayRead(da,localU,&uarray); 
   DMDAVecGetArray(da,F,&f); 
   DMDAVecGetArray(da,Udot,&udot); 
  /* Get local grid boundaries */
   DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 
  /* Compute function over the locally owned part of the grid */

  for (k=zs; k<zs+zm; k++) 
  for (j=ys; j<ys+ym; j++) 
    for (i=xs; i<xs+xm; i++) {
	    
	gk = k - zs; gj = j - ys; gi = i - xs;

	sx 		= diffusion /(DX * DX); 	
	sy 		= diffusion /(DY * DY);  
	sz 		= diffusion / (DZ * DZ);

if(user->geometry[gk][gj][gi][0]==2){ // SAN
	sx 		= (0.025 + user->dist[gk][gj][gi] ) * diffusion /(DX * DX); 	
	sy 		= (0.025 + user->dist[gk][gj][gi] ) * diffusion /(DY * DY);  
	sz 		= (0.025 + user->dist[gk][gj][gi] ) * diffusion / (DZ * DZ);

} else if(user->geometry[gk][gj][gi][0]==1){ // atrial
	sx 		= diffusion /(DX * DX); 	
	sy 		= diffusion /(DY * DY);  
	sz 		= diffusion /(DZ * DZ);
}


      /* Boundary conditions */
      if (i == 0 || j == 0 || i == usr_MX-1 || j == usr_MY-1 || k == 0 || k == usr_MZ-1 ) {
//	      ux = uy = ddx = ddy = 0.0; // is this true in uniform homogeneous atrial tissue?
/* Neumann BC */
// 6 surfaces
          	        if (i == 0 			&& j > 0 	&& j < usr_MY-1 && k > 0 && k < usr_MZ-1)    f[k][j][i] = -uarray[k][j][i] + uarray[k][j][i+1]; 
		else if (i == usr_MX-1 	&& j > 0 	&& j < usr_MY-1 && k > 0 && k < usr_MZ-1)    f[k][j][i] = -uarray[k][j][i] + uarray[k][j][i-1];  
		else if (j ==  0 		&& i > 0 	&& i < usr_MX-1 && k > 0 && k < usr_MZ-1)     f[k][j][i] = -uarray[k][j][i] + uarray[k][j+1][i]; 
		else if (j == usr_MY-1 	&& i > 0 	&& i < usr_MX-1 && k > 0 && k < usr_MZ-1)     f[k][j][i] = -uarray[k][j][i] + uarray[k][j-1][i]; 
		else if (k ==  0 		&& i > 0 	&& i < usr_MX-1 && j > 0  && j < usr_MY-1)      f[k][j][i] = -uarray[k][j][i] + uarray[k+1][j][i]; 		
		else if (k == usr_MZ-1 	&& i > 0 	&& i < usr_MX-1 && j > 0 && j < usr_MY-1)       f[k][j][i] = -uarray[k][j][i] + uarray[k-1][j][i]; 				
// 12 edges.                
else if (i == 0               && j == 0 		    && (k>0 && k<usr_MZ-1) )  	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i+1]) * sx + (-uarray[k][j][i] + uarray[k][j+1][i]) * sy ;
else if (i == 0               && j == usr_MY-1  && (k>0 && k<usr_MZ-1)) 	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i+1]) * sx + (-uarray[k][j][i] + uarray[k][j-1][i]) * sy ;
else if (i == usr_MX-1 && j == 0 		    && (k>0 && k<usr_MZ-1) )  	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i-1]) * sx + (-uarray[k][j][i] + uarray[k][j+1][i]) * sy ;
else if (i == usr_MX-1 && j == usr_MY-1  && (k>0 && k<usr_MZ-1)) 	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i-1]) * sx + (-uarray[k][j][i] + uarray[k][j-1][i]) * sy ;
else if (i == 0               && k == 0  		    && (j>0 && j<usr_MY-1))  f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i+1])*sx + (-uarray[k][j][i] + uarray[k+1][j][i]) * sz ;
else if (i == 0               && k == usr_MZ-1  && (j>0 && j<usr_MY-1))   	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i+1])*sx + (-uarray[k][j][i] + uarray[k-1][j][i]) * sz ;
else if (i == usr_MX-1 && k == 0  		    && (j>0 && j<usr_MY-1))  f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i-1])*sx + (-uarray[k][j][i] + uarray[k+1][j][i]) * sz ;
else if (i == usr_MX-1 && k == usr_MZ-1  && (j>0 && j<usr_MY-1))   	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j][i-1])*sx + (-uarray[k][j][i] + uarray[k-1][j][i]) * sz ;
else if (j == 0               && k == 0 		    && (i>0 && i<usr_MX-1) ) f[k][j][i] = (-uarray[k][j][i] + uarray[k][j+1][i])*sy + (-uarray[k][j][i] + uarray[k+1][j][i]) * sz ;
else if (j == 0               && k == usr_MZ-1  && (i>0 && i<usr_MX-1) )  	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j+1][i])*sy + (-uarray[k][j][i] + uarray[k-1][j][i]) * sz ;                
else if (j == usr_MY-1 && k == 0 		    && (i>0 && i<usr_MX-1) ) f[k][j][i] = (-uarray[k][j][i] + uarray[k][j-1][i])*sy + (-uarray[k][j][i] + uarray[k+1][j][i]) * sz ;
else if (j == usr_MY-1 && k == usr_MZ-1  && (i>0 && i<usr_MX-1) )  	f[k][j][i] = (-uarray[k][j][i] + uarray[k][j-1][i])*sy + (-uarray[k][j][i] + uarray[k-1][j][i]) * sz ;
// 8 points
else if (i == 0               && j == 0               && k == 0)   			f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i+1])*sx+(-uarray[k][j][i]+uarray[k][j+1][i])*sy+(-uarray[k][j][i]+uarray[k+1][j][i])*sz; 
else if (i == 0               && j == 0               && k == usr_MZ-1)   f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i+1])*sx+(-uarray[k][j][i]+uarray[k][j+1][i])*sy+(-uarray[k][j][i]+uarray[k-1][j][i])*sz;                 
else if (i == usr_MX-1 && j == 0               && k == 0)   		f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i-1])*sx+(-uarray[k][j][i]+uarray[k][j+1][i])*sy+(-uarray[k][j][i]+uarray[k+1][j][i])*sz;
else if (i == usr_MX-1 && j == 0               && k == usr_MZ-1) 	f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i-1])*sx+(-uarray[k][j][i]+uarray[k][j+1][i])*sy+(-uarray[k][j][i]+uarray[k-1][j][i])*sz;
else if (i == 0               && j == usr_MY-1 && k == 0)   		f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i+1])*sx+(-uarray[k][j][i]+uarray[k][j-1][i])*sy+(-uarray[k][j][i]+uarray[k+1][j][i])*sz;
else if (i == 0               && j == usr_MY-1 && k == usr_MZ-1)   f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i+1])*sx+(-uarray[k][j][i]+uarray[k][j-1][i])*sy+(-uarray[k][j][i]+uarray[k-1][j][i])*sz;                
else if (i == usr_MX-1 && j ==  usr_MY-1&& k == 0)   		f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i-1])*sx+(-uarray[k][j][i]+uarray[k][j-1][i])*sy+(-uarray[k][j][i]+uarray[k+1][j][i])*sz;
else if (i == usr_MX-1 && j ==  usr_MY-1&& k == usr_MZ-1) 	f[k][j][i]=(-uarray[k][j][i]+uarray[k][j][i-1])*sx+(-uarray[k][j][i]+uarray[k][j-1][i])*sy+(-uarray[k][j][i]+uarray[k-1][j][i])*sz;                                 
      } else { /* Interior */
        u = uarray[k][j][i];
        /* 7 point stencil */
        uxx  	= ( -2.0 * u +  uarray[k][j][i-1]  + uarray[k][j][i+1]);
        if(user->geometry[gk][gj][gi][2]==0){         uxx  	= ( -2.0 * u +  uarray[k][j][i-1]  + uarray[k][j][i-1]); }
        if(user->geometry[gk][gj][gi][4]==0){        uxx  	= ( -2.0 * u +  uarray[k][j][i+1]  + uarray[k][j][i+1]); }        
        if(user->geometry[gk][gj][gi][2]==0 && user->geometry[gk][gj][gi][4]==0) {uxx = 0.0; 			 }
        
        uyy 	= ( -2.0 * u +  uarray[k][j-1][i]  + uarray[k][j+1][i]);
        if(user->geometry[gk][gj][gi][1]==0){ uyy	= ( -2.0 * u +  uarray[k][j-1][i]  + uarray[k][j-1][i]);		 }
        if(user->geometry[gk][gj][gi][3]==0){ uyy	= ( -2.0 * u +  uarray[k][j+1][i]  + uarray[k][j+1][i]);	 }
        if(user->geometry[gk][gj][gi][1]==0 && user->geometry[gk][gj][gi][3]==0) {uyy = 0.0; 			 }
        
        uzz 		= (  -2.0 * u +  uarray[k-1][j][i]  + uarray[k+1][j][i]);
	if(user->geometry[gk][gj][gi][5]==0) {        uzz 		= (  -2.0 * u +  uarray[k-1][j][i]  + uarray[k-1][j][i]); }
	if(user->geometry[gk][gj][gi][6]==0) {        uzz 		= (  -2.0 * u +  uarray[k+1][j][i]  + uarray[k+1][j][i]); }	
	if(user->geometry[gk][gj][gi][5]==0 && user->geometry[gk][gj][gi][6]==0) uzz = 0.0;
	
        // simple
	if(user->geometry[gk][gj][gi][0]==0){ sx = sy = sz = uxx = uyy = uzz = uxy = 0.0; }
        
        // then set up the residual, f[j][i]. This f is given to the Newton's iterations.
        f[k][j][i] 	= udot[k][j][i] - (uxx * sx + uyy * sy + uzz * sz );
      }
    } // end of the k, j, i loops.

  /* Restore vectors */
   DMDAVecRestoreArrayRead(da,localU,&uarray); 
   DMDAVecRestoreArray(da,F,&f); 
   DMDAVecRestoreArray(da,Udot,&udot); 
   DMRestoreLocalVector(da,&localU); 
  PetscFunctionReturn(0);
}

