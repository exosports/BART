/*
 * geometry.c   - Functions to establish a system geometry. Component
 *                of the transit program
 *
 * Copyright (C) 2003-2006 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of version 2 of the GNU General 
 * Public License as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include <transit.h>

/*\fcnfh
  Takes hinted geometry structure and fill it up with calculated values.

  @returns 0 on success
*/
int
setgeomhint(struct transit *tr)
{
  static struct geometry st_sg;
  memset(&st_sg,0,sizeof(struct geometry));

  struct geometry *sg=tr->ds.sg=&st_sg;
  struct geometry *hg=&tr->ds.th->sg;

  sg->transpplanet=hg->transpplanet;

  //If factor values are correctly hinted then use them, otherwise
  //assume hard coded values. If you change the following lines, don't
  //forget to change the help page in transit.h that the semimajor axis
  //is in AU; timing is in hours from mid eclipse; star's mass and
  //radius are in solar units.
  sg->smaxisfct  =hg->smaxisfct>0  ?hg->smaxisfct  :AU;
  sg->timefct    =hg->timefct>0    ?hg->timefct    :HOUR;
  sg->eccfct     =hg->eccfct>0     ?hg->eccfct     :1;
  sg->inclfct    =hg->inclfct>0    ?hg->inclfct    :DEGREES;
  sg->aperfct    =hg->aperfct>0    ?hg->aperfct    :DEGREES;
  sg->lnode      =hg->lnodefct>0   ?hg->lnodefct   :DEGREES;

  sg->starmassfct=hg->starmassfct>0?hg->starmassfct:SUNMASS;
  sg->starradfct =hg->starradfct>0 ?hg->starradfct :SUNRADIUS;


  //If values are not hinted, then use hard coded values: 
  sg->smaxis  =hg->smaxis>0  ?hg->smaxis  :AU/sg->smaxisfct;
  sg->time    =hg->time>0    ?hg->time    :0;
  sg->ecc     =hg->ecc>0     ?hg->ecc     :0;
  sg->incl    =hg->incl>0    ?hg->incl    :0;
  sg->aper    =hg->aper>0    ?hg->aper    :0;
  sg->lnode   =hg->lnode>0   ?hg->lnode   :0;

  sg->starmass=hg->starmass>0?hg->starmass:1.101*SUNMASS/sg->starmass;
  sg->starrad =hg->starrad>0 ?hg->starrad :1.125*SUNRADIUS/sg->starradfct;

  //set progressindicator and return
  tr->pi|=TRPI_GEOMETRYHINT;
  return 0;
}


/* \fcnfh
   Set geometry variables

   @returns 0 on success
*/
int
setgeom(struct geometry *sg,	/* geometry structure */
	double time,		/* time for which to set planet
				   position, use HUGE_VAL for use the
				   default value stored in the geometry
				   structure */
	long *flags)		/* Progress indicator flag */
{
  transitcheckcalled(*flags,"setgeom",1,
		     "setgeomhint",TRPI_GEOMETRYHINT
		     );

  /* TD: include argument of pericenter and longitud of the node in the
     computations, right now both values are used as zero. */
  //Auxiliary variables in cgs
  double smaxis=sg->smaxis*sg->smaxisfct;
  double ecc=sg->ecc*sg->eccfct;
  double incl=sg->incl*sg->inclfct;
  //  double aper=sg->aper*sg->aperfct;
  //  double lnode=sg->lnode*sg->lnodefct;
  double t=(time<HUGE_VAL?time:sg->time)*sg->timefct;
  double mass=sg->starmass*sg->starmassfct;

  //set mean angular velocity, mean true anomaly, and others variables
  double prec2=1e-12;
  double n=sqrt(GGRAV*mass/smaxis/smaxis/smaxis);
  double Ea=HUGE_VAL,E=n*time;
  while((E-Ea)*(E-Ea)>prec2){
    Ea=E;
    E=n*t+ecc*sin(E);
  }
  double M=smaxis*(1-ecc*ecc);
  double Delta=smaxis*(1-ecc*cos(E));
  double cosv=(M-Delta)/ecc/Delta;
  double sini=sin(incl);
  double cosi=cos(incl);

  //set X and Y.
  sg->x=Delta*sqrt(cosi*cosi-cosv*cosv);
  sg->y=Delta*sini;
  
  //set progress indicator and return.
  *flags|=TRPI_GEOMETRY;
  return 0;
}


/* \fcnfh
   Any normalized variation of the star, probably limb darkening

   @return normalized value;
*/
inline PREC_RES
starvariation(double x,
	      double y,
	      double radius)
{
  if(x*x+y*y>radius*radius)
    return 0;

  return 1;
}
