/*
 * geometry.c
 * geometry.txc - Functions to establish a system geometry. Component
 *                of the transit program
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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
setgeom(struct transit *tr)
{
  static struct geometry st_sg;
  memset(&st_sg,0,sizeof(struct geometry));

  struct geometry *sg=tr->ds.sg=&st_sg;
  struct geometry *hg=&tr->ds.th->sg;


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

  sg->starmass=hg->starmass>0?hg->starmass:SUNMASS/sg->starmass;
  sg->starrad =hg->starrad>0 ?hg->starrad :SUNRADIUS/sg->starradfct;


  tr->pi|=TRPI_GEOMETRY;
  return 0;
}
