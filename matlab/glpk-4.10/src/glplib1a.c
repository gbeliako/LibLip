/* glplib1a.c (platform-independent ISO C version) */

/*----------------------------------------------------------------------
-- This code is part of GNU Linear Programming Kit (GLPK).
--
-- Copyright (C) 2000, 01, 02, 03, 04, 05, 06 Andrew Makhorin,
-- Department for Applied Informatics, Moscow Aviation Institute,
-- Moscow, Russia. All rights reserved. E-mail: <mao@mai2.rcnet.ru>.
--
-- GLPK is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2, or (at your option)
-- any later version.
--
-- GLPK is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
-- or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
-- License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with GLPK; see the file COPYING. If not, write to the Free
-- Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
-- 02110-1301, USA.
----------------------------------------------------------------------*/

#include <stddef.h>
#include "glplib.h"

static void *pointer = NULL;
/* some "secret" place to store a pointer */

/*----------------------------------------------------------------------
-- lib_set_ptr - store a pointer.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- void lib_set_ptr(void *ptr);
--
-- *Description*
--
-- The routine lib_set_ptr stores a pointer specified by the parameter
-- ptr in some "secret" place. */

void lib_set_ptr(void *ptr)
{     pointer = ptr;
      return;
}

/*----------------------------------------------------------------------
-- lib_get_ptr - retrieve a pointer.
--
-- *Synopsis*
--
-- #include "glplib.h"
-- void *lib_get_ptr(void);
--
-- *Returns*
--
-- The routine lib_get_ptr returns a pointer previously stored by the
-- routine lib_set_ptr. If the latter has not been called yet, NULL is
-- returned. */

void *lib_get_ptr(void)
{     void *ptr;
      ptr = pointer;
      return ptr;
}

/* eof */
