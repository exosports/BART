/*********************************************************************
 * 
 * dbread_debug.c - Driver to read a dataset which is intended to debug
 *                  lineread.
 *
 * Copyright (C) 2006 Patricio Rojo
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 *
 **********************************************************************/

#include <lineread.h>


static const driver_func pdriverf_debug = {
  "DEBUGGING driver",
  NULL,
  &debug_find,
  &debug_open,
  &debug_close,
  &debug_info,
  &debug_partition
};

const driver_func *driverf_debug = &pdriverf_debug;
