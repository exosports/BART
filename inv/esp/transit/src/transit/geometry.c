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
  struct geometry *sg=tr->ds.sg=&st_sg;
  struct geometry *hg=&tr->ds.th->sg;

  sg->smaxis=hg->smaxis;

  //If values are not hinted, then use hard coded values: semimajor axis
  //of 1AU. Timing is mid eclipse. star's mass and radius are solar
  sg->smaxis=hg->smaxis>0?hg->smaxis:1*AU;
  sg->time=hg->time>0?hg->time:0;
  sg->starmass=hg->starmass>0?hg->starmass:1*SUNMASS;
  sg->starrad=hg->starrad>0?hg->starrad:1*SUNRADIUS;

  //If factor values are correctly hinted then use them, otherwise just
  //assume that values are in cgs.
  sg->smaxisfct=hg->smaxisfct>0?hg->smaxisfct:1;
  sg->timefct=hg->timefct>0?hg->timefct:1;
  sg->starmassfct=hg->starmassfct>0?hg->starmassfct:1;
  sg->starradfct=hg->starradfct>0?hg->starradfct:1;

  return 0;
}
