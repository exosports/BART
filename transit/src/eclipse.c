
#include <transit.h>
#include <extraext.h>

/* initial version March 2nd, 2014 Jasmina Blecic                           */



/* ##########################
    CALCULATES OPTICAL DEPTH
   ########################## */

/* \fcnfh
   Computes optical depth for eclipse geometry for one ray, one wn, 
   between a certain layer in the atmosphere up to the top layer. 
   Returns: Optical depth divided by rad.fct:  \frac{tau}{units_{rad}}      */

static PREC_RES
totaltau_eclipse( PREC_RES *rad,   /* Equispaced radius array               */
                                   /* Current layer to outmost layer        */
                  PREC_RES *ex,    /* Extinction[rad], all radii, one wn    */
                                   /* Often used as ex+ri than &ex[ri]      */
                    long nrad)    /* Number of radii                       */
                                   /* Current layer to outmost layer        */
{
  PREC_RES res;         /* Optical depth divided by units of radius         */
  PREC_RES x3[3],r3[3]; /* Interpolation variables                          */

  /* Distance between two radii. Radius array needs to be equispaced.       */
  const PREC_RES dr=rad[1]-rad[0];

  /* Distance along the path                                                */
  PREC_RES s[nrad];

  /* Returns 0 if at outmost layer. No distance travelled.                  */
  if(nrad == 1)  
    return 0.;

  /* Providing three necessary points for spline integration:               */
  const PREC_RES tmpex=*ex;
  const PREC_RES tmprad=*rad;
  if(nrad==2) *ex=interp_parab(rad-1,ex-1, rad[0]);
  else *ex=interp_parab(rad,ex, rad[0]);

  if(nrad==2){
    x3[0]=ex[0];
    x3[2]=ex[1];
    x3[1]=(ex[1]+ex[0])/2.0;
    r3[0]=rad[0];
    r3[2]=rad[1];
    r3[1]=(rad[0]+rad[1])/2.0;
    *rad=tmprad;
    *ex=tmpex;
    rad=r3;
    ex=x3;
    nrad++;
  }

  /* Distance along the path:                                               */
  s[0] = 0.;
  for(int i=1; i< nrad; i++){
    s[i]= s[i-1] + dr;          
  }

  /* Integrate extinction along the path:                                   */
  /* Use spline if GSL is available along with at least 3 points:           */
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spl = gsl_interp_alloc(gsl_interp_cspline, nrad);
  gsl_interp_init(spl, s, ex, nrad);
  res = gsl_interp_eval_integ(spl, s, ex, 0, s[nrad-1], acc);    
  gsl_interp_free(spl);
  gsl_interp_accel_free(acc);
#else
#error non equispaced integration is not implemented without GSL
#endif /* _USE_GSL */

  /* Optical depth divided by units of radius:                              */
  return res;
}


/* ################################### 
    FILLS OUT OPTICAL DEPTH 2D ARRAY
   ################################### */



#define CIA_DOFLOAT  2
#define CIA_RADFIRST 1

/* \fcnfh
   Calculates optical depth as a function of radii and wavenumber.
   Return: 0 on success                                                     */

int
tau_eclipse(struct transit *tr){  /* Transit structure             */
  static struct optdepth tau;              /* Def optical depth structure   */
  tr->ds.tau = &tau;                       /* Called from transit structure */
  prop_samp *rad = &tr->rads;              /* Radius sample                 */
  prop_samp *wn  = &tr->wns;               /* Wavenumber sample             */
  struct extinction *ex = tr->ds.ex;       /* Def extinction structure      */
  PREC_RES **e = ex->e[tr->tauiso];        /* 2D extinction array           */
                                           /* One isotope, [rad][wn]        */
  PREC_RES (*fcn)() = tr->ecl->tauEclipse; /* Defines function fcn          */
                                           /* Calls totaltau_eclipse        */
                                           /* See structure_tr.h line 91    */

  long wi, ri; /* Indices for wavenumber, and radius                       */
  int rn;      /* Returns error code from computeextradiaus=>exstradis     */

  PREC_RES *r  = rad->v;         /* Values of radii                         */
  PREC_RES *tau_wn;              /* Array of optical depths [rad], one wn   */
  PREC_ATM *temp = tr->atm.t,    /* Temp array from atmospheric file        */
           tfct  = tr->atm.tfct; /* Temp units from atmospheric file        */

  /* Number of elements:                                                   */
  long int wnn = wn->n;        /* Wavenumbers                              */
  long int rnn = rad->n;       /* Radii                                    */
  double wfct  = wn->fct;      /* Wavenumber units factor to cgs           */

  /* Checks weather extinction exists:                                     */
  transitcheckcalled(tr->pi, "tau", 1, "extwn", TRPI_EXTWN);


  /* Pato Rojo explanation: 
  Consider:
  TRU_*BITS
  Those mark the bits that are carrying flags for each aspect. 
  i.e. TRU_EXTBITS is 1 for all the flag positions that are relevant
  for the extinction computation (only TRU_OUTTAU in that case); 
  TRU_ATMBITS is 1 for all the flag positions that are relevant for 
  the atmospheric computation; and so  on...therefore the line 64 
  passes all the TAU relevant flags from the user hint to 
  the main structure.  In particular, it only passes whether TRU_OUTTAU 
  was 1 or 0 as specified by the user. */
  transitacceptflag(tr->fl, tr->ds.th->fl, TRU_TAUBITS);


  /* Sets maximum optical depth:                                            */
  struct transithint *trh = tr->ds.th;
  tau.toomuch = 10;            /* Default value if not set in conf file     */
  if(tr->ds.th->toomuch > 0)   
    tau.toomuch = trh->toomuch; 

  /* Allocates array to store radii indices where tau reaches toomuch:      */
  tau.last = (long      *)calloc(wnn,       sizeof(long));
  /* Allocates 2D array [rad][wn]:                                          */
  tau.t    = (PREC_RES **)calloc(wnn,        sizeof(PREC_RES *));
  tau.t[0] = (PREC_RES  *)calloc(wnn*rad->n, sizeof(PREC_RES  ));

  /* Connects tau.t with correct wn and rad:                                */
  for(wi=1; wi < wnn; wi++)
    tau.t[wi] = tau.t[0] + wi*rad->n;

  /* Set cloud structure:                                                   */
  /* Linear for rini to rfin, grey opacity from 0 to maxe.                  */
  static struct extcloud cl;
  cl.maxe = tr->ds.th->cl.maxe; /* Maximum opacity                          */
  cl.rini = tr->ds.th->cl.rini; /* Top layer radius                         */
  cl.rfin = tr->ds.th->cl.rfin; /* Radius of maxe                           */
  if(tr->ds.th->cl.rfct==0)
    cl.rfct = rad->fct;
  else
    cl.rfct = tr->ds.th->cl.rfct;
  tr->ds.cl = &cl;

  /* Tests if ex is calculated for particular radius, if true/False:        */
  _Bool *comp = ex->computed;

  /* Restors extinction from the savefile if given:                         */
  if(tr->save.ext)
    restfile_extinct(tr->save.ext, e, comp, rnn, wnn);


  /* Computes extinction at the outermost layer:                            */
  /* Pato Rojo explanation:                                                 */
  /* If outmost layer not computed r[lastr] will segfault otherwise.        */
  /* Last computed layer in that case is unnecessary.                       */
  /* Compute extinction when radius is below the last computed layer:       */
  if(!comp[rnn-1]){
    transitprint(1, verblevel,
                 "Computing extinction in the outtermost layer.\n");
    if((rn=computeextradius(rnn-1, tr->atm.t[rnn-1]*tr->atm.tfct, ex))!=0)
      transiterror(TERR_CRITICAL,
                   "computeextradius() returned error code %i.\n", rn);
  }

  /* Gets a copy of the radius units factor:                                */
  double rfct = rad->fct;

  /* Requests at least four radii to calculate a spline interpolation:      */
  if(rnn < 4)
    transiterror(TERR_SERIOUS, "tau(): At least four radius points are "
                               "required! (three for spline and one for the "
                               "analytical part)");

  transitprint(1, verblevel,
                  "Calculating optical depth at various radius ...\n");

  /* Note that it works only for one isotope:                               */
  if(ex->periso)
    transitprint(2, verblevel,
                 "Computing only for isotope '%s', others were ignored.\n",
                 tr->ds.iso->isof[tr->tauiso].n);

  PREC_RES er[rnn];        /* Array of extinction per radius                */
  int lastr = rnn-1;       /* Radius index of last computed extinction      */
                           /* Starts from the outmost layer                 */

  /* Tenth of wavenumber sample size,used for progress printing:            */
  int wnextout = (long)(wnn/10.0); 


  /* Extinction from scattering and clouds:                                 */
  double e_s[rnn], e_c[rnn];
  /* Extinction from CIA:                                                   */
  PREC_CIA **e_cia = tr->ds.cia->e;
  /* Extinction from scattering:                                            */
  /* This is a hook for future handling of scattering.                      */
  struct extscat *sc = tr->ds.sc;



  /* Flow of this part:
     For each wavenumber of the selected scale, the optical depth along 
     the path with a given radius is computed. 
     The computation starts with an radius r[rnn], equal to the radius 
     of the outmost layer and descends into the atmosphere until the 
     optical depth rises above a certain user-defined limit (--toomuch). 
     Each time this procedure reaches a radius whose extinction has not 
     been computed, it calls the extinction computing routine, 
     to compute the extinction over the whole wavenumber range .            */
  

  /* For each wavenumber:                                                   */
  for(wi=0; wi < wnn; wi++){
    /* Pointing to the declared tau array. All radii for one wn.            */
    tau_wn = tau.t[wi];

    /* Print output every 10\% that is ready:                               */
    if(wi > wnextout){
      transitprint(2, verblevel, "%i%%\n", (int)(100*(float)wi/wnn+0.5));
      wnextout += (long)(wnn/10.0);
    }

    /* Calculate extinction from scattering, clouds, and CIA.               */
    /* At each level just for the wn of interest:                           */
    computeextscat(e_s,  rnn, sc, rad->v, rad->fct, temp, tfct, wn->v[wi]*wfct);
    computeextcloud(e_c, rnn, &cl, rad, temp, tfct, wn->v[wi]*wfct);

    /* Calculate extinction for all radii if extinction exists:             */
    for(ri=0; ri < rnn; ri++)
      er[ri] = e[ri][wi] + e_s[ri] + e_c[ri] + e_cia[wi][ri];

    /* For each radii:                                                      */
    for(ri=rnn-1; ri > -1; ri--){
      transitprint(10, verblevel, "Radius r[%li]=%9.4g\n", ri, r[ri]);

      /* Computes extinction only if radius is smaller then 
         the radius of the last calculated extinction.                      */
      if(r[ri] < r[lastr]){    
        transitprint(3, verblevel, "Last Tau (r=%9.4g, wn=%9.4g): %10.4g.\n",
                                     r[ri-1], wn->v[wi], tau_wn[ri-1]);

        /* If extinction was not computed, compute it using computextradius. 
           The function calls extradius for one radius, but for all wn.
           Then, the code updates extinction for the particular wn. 
           Computation starts from the layer below outermost, since the 
           outmost layer is computed above.                                 */ 
                                          
        do{
          if(!comp[--lastr]){
            /* Compute a new extinction at given radius printing error if
               something happen:                                            */
            transitprint(2, verblevel, "Radius %i: %.9g cm ... ",
                                       lastr+1, r[lastr]);
            if((rn=computeextradius(lastr, temp[lastr]*tfct, ex))!=0)
              transiterror(TERR_CRITICAL,
                           "computeextradius() return error code %i while "
                           "computing radius #%i: %g\n", rn, r[lastr]*rfct);
            /* Update the value of the extinction at the right place:       */
            er[lastr] = e[lastr][wi] + e_s[lastr] +
                        e_c[lastr] + e_cia[wi][lastr];
          }
        }while(r[ri] < r[lastr]);
      }

      /* Fills out tau_wn[ri] array until tau reaches toomuch
         first tau[0], starts from the outmost layer.                         
         Also fills out tau.last array, which gives radii
         where tau reaches toomuch for each wn.                             */ 
      if( (tau_wn[rnn-ri-1] = rfct * fcn(r+ri, er+ri, rnn-ri)) > tau.toomuch){
        tau.last[wi] = rnn-ri-1;

        if (ri < 3){
          transitprint(1, verblevel,
                       "WARNING: At wavenumber %g (cm-1), the critical TAU "
                       "value (%g) was exceeded with tau=%g at the radius "
                       "level %li (%g km), this should have "
                       "happened in a deeper layer (check IP sampling or ATM "
                       "file).\n", wn->v[wi], tau.toomuch, tau_wn[ri],
                       ri, r[ri]*rfct/1e5);
        }
      /* Exit impact-parameter loop if it reached toomuch:                  */
        break;
      }

      /* Sets tau of the outermost layer to zero:                           */
      tau_wn[0] = 0;   

      transitDEBUG(22, verblevel,
                   "Tau(lambda %li=%9.07g, r=%9.4g) : %g  (toomuch: %g)\n",
                   wi, wn->v[wi], r[ri], tau_wn[ri], tau.toomuch);
    }

    if(ri==rnn){
      transitprint(1, verblevel,
                   "WARNING: At wavenumber %g cm-1, the bottom of the "
                   "atmosphere was reached before obtaining the critical TAU "
                   "value of %g.\nMaximum TAU reached: %g.\n",
                   wn->v[wi], tau.toomuch, tau_wn[ri]);
      tau.last[wi] = ri-1;
    }

  }

  transitprint(1, verblevel, " Done.\nOptical depth calculated up to "
                             "tau=%g.\n", tr->ds.tau->toomuch);

  /* Print detailed output if requested:                                    */
  if(tr->ds.det->tau.n)
    detailout(&tr->wns, &tr->rads,  &tr->ds.det->tau, tau.t, 0);
  if(tr->ds.det->ext.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->ext, e, CIA_RADFIRST);
  if(tr->ds.det->cia.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->cia, (double **)e_cia,
              CIA_DOFLOAT);

  if(tr->save.ext)
    savefile_extinct(tr->save.ext, e, comp, rnn, wnn);


  /* Print lowest radius before optical depth gets too big:                 */
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch, tr->ds.tau, &tr->wns, &tr->rads);

  /* Free memory that is no longer needed:                                  */
  //freemem_lineinfotrans(tr->ds.li, &tr->pi);
  freemem_localextinction();

  /* Set progress indicator and output tau if requested, 
     otherwise return success:                                              */
  tr->pi |= TRPI_TAU;
  if(tr->fl & TRU_OUTTAU)
    printtau(tr);
  return 0;
}

/* ################################################# 
    CALCULATES EMERGENT INTENSITY FOR ONE WAVENUMBER
   ################################################# */


/* \fcnfh
   Calculates emergent intensity.
   Return: emergent intensity for one wavenumber                            */

static PREC_RES
eclipse_intens(struct transit *tr,  /* Transit structure                   */
            PREC_RES *tau,          /* Optical depth array tau.c==>tau     */
            PREC_RES w,             /* Current wavenumber value            */
            long last,              /* Index where tau = toomuch           */
            double toomuch,         /* Maximum optical depth calculated    */
            prop_samp *rad){        /* Radii array                         */


  /* General variables:                                                     */
  PREC_RES res;                  /* Result                                  */
  PREC_ATM *temp = tr->atm.t;    /* Temperatures                            */

  /* Takes sampling properties for wavenumber from tr:                      */
  prop_samp *wn = &tr->wns; 
  /* Wavenumber units factor to cgs:                                        */
  double wfct  = wn->fct; 

  /* Radius parameter variables:                                            */
  long rnn  = rad->n;
  long i;

  /* Blackbody variables:                                                   */
  PREC_RES B[rnn];

  /* Integration parts:                                                     */
  PREC_RES tauInteg[rnn],  /* Integrand function                            */
           tauIV[rnn];     /* Tau integration variable                      */

  /* Integrate for each of the planet's layer starting from the           
     outermost until the closest layer. 
     The order is opposite since tau starts from the top and 
     radius array starts from the bottom.                                   */

  /* Plank function for wavenumbers:
     B_wn =\frac{2.*PlankConst*wn^3*c^2}
                {exp^(\frac{PlankConst*wn*c}{Kb*T})-1.}                     
     Units for Plunk Function are: [erg/s/sr/cm].                           */
  for(i=0;i <= last; i++){   
    tauIV[i] = tau[i];  
    B[i] =  (2. * H * w * w * w * wfct * wfct * wfct * LS * LS) 
          / (exp(H * w * wfct * LS / (KB * temp[rnn-1-i])) - 1.);         
    tauInteg[i] = B[i] * exp(-tau[i]);
  }


  /* Added 0 at the end when tau reach toomuch, so the spline looks nice    */
  /* Add all other layers to be 0.                                          */
  for(; i<rnn; i++){
    tauInteg[i] = 0;
    /* Geometric progression is used to provide enough elements 
       for integral to work. It does not change the final outcome/result.   */
    tauIV[i]= tauIV[i-1] + 1; 
   }


  /* Adding additional 0 layer, plus the last represent number of elements 
     is -1, so we need to add one more. 2 in total.                         */
  last+=2;  

  /* If atmosphere is transparent, and at last level tau has not reached 
     tau.toomuch, last is set to max number of layers (rnn, instead of rnn-1
     because we added 2 on the previous step). The code requests never
     to go over it.                                                         */
  if(last > rnn)    
    last= rnn;

  /* Checks if we have enough radii to do spline, at least 3:               */
  if(last < 3)
    transiterror(TERR_CRITICAL,
                 "Condition failed, less than 3 items (only %i) for radial "
                 "integration.\n", last);

  /* Integrate along tau up to tau = toomuch:                               */
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spl       = gsl_interp_alloc(gsl_interp_cspline, last);
  gsl_interp_init(spl, tauIV, tauInteg, last);
  res = gsl_interp_eval_integ(spl, tauIV, tauInteg,
                               tauIV[0], tauIV[last-1], acc);
  gsl_interp_free(spl);
  gsl_interp_accel_free (acc);
#else
# error computation of modulation() without GSL is not implemented
#endif

  return res;
}



/* ############################################################# 
    CALCULATES EMERGENT INTENSITY FOR THE WHOLE WAVENUMBER RANGE
   ############################################################# */


/* \fcnfh
   Calculates the emergent intensity for the whole range of wavenumbers
   Returns: emergent intensity for the whole wavenumber range               */

int
emergent_intens(struct transit *tr){  /* Transit structure                 */

  static struct outputray st_out;    /* Defines output structure          */
  tr->ds.out = &st_out;


  /* Initial variables:                                                     */
  long w;
  prop_samp *rad = &tr->rads;          /* Radius array pointer              */
  prop_samp *wn = &tr->wns;            /* Wavenumber array pointer          */
  eclipse_ray_solution *ecl = tr->ecl; /* Eclipse ray solution pointer      */


  /* Allocate the intensity array:                                          */
  PREC_RES *out = st_out.o = (PREC_RES *)calloc(wn->n, sizeof(PREC_RES));
  struct optdepth *tau = tr->ds.tau;


  /* Integrate for each wavelength:                                         */
  transitprint(1, verblevel, "\nIntegrating over wavelength.\n");

  /* Printing process variable:                                             */
  int nextw = wn->n/10;


  /* Calculate the intensity integral at each wavenumber:                   */
  for(w=0; w<wn->n; w++){
    out[w] = ecl->eclIntenWn(tr, tau->t[w], wn->v[w], tau->last[w], 
                             tau->toomuch, rad);

    /* Print to screen the progress status:                                 */
    if(w==nextw){
      nextw += wn->n/10;
      transitprint(2, verblevel, "%i%% ", (10*(int)(10*w/wn->n+0.9999999999)));
    }
  }
  transitprint(1, verblevel, "\nDone.\n");

  /* Free no longer needed memory:                                          */
  freemem_idexrefrac(tr->ds.ir, &tr->pi);
  freemem_extinction(tr->ds.ex, &tr->pi);
  freemem_tau(tr->ds.tau,       &tr->pi);

  /* Set progress indicator, and print output:                              */
  tr->pi &= TRPI_MODULATION;
  printecl(tr);  
  return 0;
}

/* #########################
    OUTPUT FILE FOR ECLIPSE
   ######################### */


/* \fcnfh
   Print (to file or stdout) the emergent intensity as function of wavenumber 
   (and wavelength)                                                          */

void
printecl(struct transit *tr){
  FILE *outf=stdout;
  struct outputray *outray=tr->ds.out;
  int rn;

  /* Opens a file:                                                          */
  if(tr->f_out&&tr->f_out[0]!='-')
    outf=fopen(tr->f_out,"w");

  transitprint(1,verblevel,
	       "\nPrinting intensity for requested conditions in '%s'\n"
	       ,tr->f_out?tr->f_out:"standard output");

  /* Prints:                                                                */
  char wlu[20],wnu[20];        /* Wavelength units name                     */
  long nsd=(long)(1e6);        /* Wavenumber units name                     */

  /* Get wavenumber units name:                                             */
  if((long)(nsd*tr->wns.fct)==nsd) strcpy(wnu,"cm");
  else if((long)(nsd*1e-1*tr->wns.fct)==nsd) strcpy(wnu,"mm");
  else if((long)(nsd*1e-4*tr->wns.fct)==nsd) strcpy(wnu,"um");
  else if((long)(nsd*1e-7*tr->wns.fct)==nsd) strcpy(wnu,"nm");
  else if((long)(nsd*1e-8*tr->wns.fct)==nsd) strcpy(wnu,"a");
  else sprintf(wnu,"%6.1g cm",1/tr->wns.fct);

  /* Get wavenumber units name:                                              */
  if((long)(nsd*tr->wavs.fct)==nsd) strcpy(wlu,"cm");
  else if((long)(nsd*1e1*tr->wavs.fct)==nsd) strcpy(wlu,"mm");
  else if((long)(nsd*1e4*tr->wavs.fct)==nsd) strcpy(wlu,"um");
  else if((long)(nsd*1e7*tr->wavs.fct)==nsd) strcpy(wlu,"nm");
  else if((long)(nsd*1e8*tr->wavs.fct)==nsd)  strcpy(wlu,"a");
  else sprintf(wlu,"%8.1g cm",tr->wavs.fct);

  
  /* Prints the header:                                                      */
  fprintf(outf,
	  "#wvl [um]%*sEmergIntens [erg/s/cm/sr]\n"
	  ,6 ," ");


  /* Prints wavenumber, wavelength, and emergent intensity:                  */
  for(rn=0;rn<tr->wns.n;rn++)
    fprintf(outf,"%-15.5g%-18.9g\n"
	    ,(1/(tr->wns.v[rn]/tr->wns.fct))*1e4,
	    outray->o[rn]);

  fclose(outf);

  return;
}


/* ############################
    CALLS ECLIPSE RAY SOLUTION
   ############################ */


  /* Executes eclipse module:                                               */
const eclipse_ray_solution eclipsepath = {
       "Eclipse Path",   /* Name of the solution                            */
       "eclipse.c",      /* Source code file name                           */
       &totaltau_eclipse,/* Per per wavenumber value computation            */
       &eclipse_intens,  /* Per wavenumber value computation                */
       };
