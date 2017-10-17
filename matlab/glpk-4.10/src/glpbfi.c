/* glpbfi.c */

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
typedef struct BFI BFI;
#define _GLPBFI_DEFINED
#include "glpbfi.h"
#include "glpinv.h"
#include "glplib.h"

struct BFI
{     int m_max;
      INV *inv;
};

BFI *bfi_create_binv(void)
{     /* create factorization of the basis matrix */
      BFI *binv;
      binv = umalloc(sizeof(BFI));
      binv->m_max = 0;
      binv->inv = NULL;
      return binv;
}

int bfi_factorize(BFI *binv, int m,
      void *info, int (*col)(void *info, int j, int rn[], double bj[]))
{     /* compute factorization of the basis matrix */
      insist(m > 0);
      if (binv->m_max < m)
      {  if (binv->inv != NULL) inv_delete(binv->inv);
         binv->m_max = m + 100;
         binv->inv = inv_create(binv->m_max, 50);
      }
      binv->inv->m = m;
      binv->inv->luf->n = m;
      return inv_decomp(binv->inv, info, col);
}

void bfi_ftran(BFI *binv, double x[], int save)
{     /* perform forward transformation (FTRAN) */
      inv_ftran(binv->inv, x, save);
      return;
}

void bfi_btran(BFI *binv, double x[])
{     /* perform backward transformation (BTRAN) */
      inv_btran(binv->inv, x);
      return;
}

int bfi_update_binv(BFI *binv, int j)
{     /* update factorization of the basis matrix */
      return inv_update(binv->inv, j);
}

void bfi_delete_binv(BFI *binv)
{     /* delete factorization of the basis matrix */
      if (binv->inv != NULL) inv_delete(binv->inv);
      ufree(binv);
      return;
}

/* eof */
