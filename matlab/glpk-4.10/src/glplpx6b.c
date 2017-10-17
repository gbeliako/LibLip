/* glplpx6b.c (interior-point solver routine) */

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

#include <string.h>
#include "glpipm.h"
#include "glplib.h"
#include "glplpx.h"

/*----------------------------------------------------------------------
-- lpx_interior - easy-to-use driver to the interior-point method.
--
-- *Synopsis*
--
-- #include "glplpx.h"
-- int lpx_interior(LPX *lp);
--
-- *Description*
--
-- The routine lpx_interior is intended to find optimal solution of an
-- LP problem, which is specified by the parameter lp.
--
-- Currently this routine implements an easy variant of the primal-dual
-- interior-point method based on Mehrotra's technique.
--
-- This routine transforms the original LP problem to an equivalent LP
-- problem in the standard formulation (all constraints are equalities,
-- all variables are non-negative), calls the routine ipm1 for solving
-- the transformed problem, and then transforms an obtained solution to
-- the solution of the original problem.
--
-- *Returns*
--
-- The routine lpx_interior returns one of the following exit codes:
--
-- LPX_E_OK       the problem has been successfully solved.
--
-- LPX_E_FAULT    the solver can't start the search because the problem
--                has no rows and/or columns.
--
-- LPX_E_NOFEAS   the problem has no feasible (primal/dual) solution.
--
-- LPX_E_NOCONV   very slow convergence or divergence.
--
-- LPX_E_ITLIM    iterations limit exceeded.
--
-- LPX_E_INSTAB   numerical instability on solving Newtonian system.
--
-- Should note that additional exit codes may appear in future versions
-- of this routine. */

struct dsa
{     /* working area used by lpx_interior routine */
      LPX *lp;
      /* original LP instance to be solved */
      int orig_m;
      /* number of rows in original LP */
      int orig_n;
      /* number of columns in original LP */
      int *ref; /* int ref[1+orig_m+orig_n]; */
      /* this array contains references from components of original LP
         to components of transformed LP */
      /*--------------------------------------------------------------*/
      /* following components define transformed LP in standard form:
         minimize    c * x + c[0]
         subject to  A * x = b,   x >= 0 */
      int m;
      /* number of rows (constraints) in transformed LP */
      int n;
      /* number of columns (variables) in transformed LP */
      int size;
      /* size of arrays ia, ja, ar (enlarged automatically) */
      int ne;
      /* number of non-zeros in transformed constraint matrix */
      int *ia; /* int ia[1+size]; */
      int *ja; /* int ja[1+size]; */
      double *ar; /* double ar[1+size]; */
      /* transformed constraint matrix in storage-by-indices format
         which has m rows and n columns */
      double *b; /* double b[1+m]; */
      /* right-hand sides; b[0] is not used */
      double *c; /* double c[1+n]; */
      /* objective coefficients; c[0] is constant term */
      /*--------------------------------------------------------------*/
      /* following arrays define solution of transformed LP computed by
         interior-point solver */
      double *x; /* double x[1+n]; */
      /* values of primal variables */
      double *y; /* double y[1+m]; */
      /* values of dual variables for equality constraints */
      double *z; /* double z[1+n]; */
      /* values of dual variables for non-negativity constraints */
};

static void calc_mn(struct dsa *dsa)
{     /* calculate the number of constraints (m) and variables (n) for
         transformed LP */
      LPX *lp = dsa->lp;
      int orig_m = dsa->orig_m;
      int orig_n = dsa->orig_n;
      int i, j, m = 0, n = 0;
      for (i = 1; i <= orig_m; i++)
      {  switch (lpx_get_row_type(lp, i))
         {  case LPX_FR:                  break;
            case LPX_LO: m++; n++;        break;
            case LPX_UP: m++; n++;        break;
            case LPX_DB: m += 2; n += 2;  break;
            case LPX_FX: m++;             break;
            default: insist(lp != lp);
         }
      }
      for (j = 1; j <= orig_n; j++)
      {  switch (lpx_get_col_type(lp, j))
         {  case LPX_FR: n += 2;          break;
            case LPX_LO: n++;             break;
            case LPX_UP: n++;             break;
            case LPX_DB: m++; n += 2;     break;
            case LPX_FX:                  break;
            default: insist(lp != lp);
         }
      }
      dsa->m = m;
      dsa->n = n;
      return;
}

static void new_coef(struct dsa *dsa, int i, int j, double aij)
{     /* add new element to the transformed constraint matrix */
      insist(1 <= i && i <= dsa->m);
      insist(1 <= j && j <= dsa->n);
      insist(aij != 0.0);
      if (dsa->ne == dsa->size)
      {  int *ia = dsa->ia;
         int *ja = dsa->ja;
         double *ar = dsa->ar;
         dsa->size += dsa->size;
         dsa->ia = ucalloc(1+dsa->size, sizeof(int));
         memcpy(&dsa->ia[1], &ia[1], dsa->ne * sizeof(int));
         ufree(ia);
         dsa->ja = ucalloc(1+dsa->size, sizeof(int));
         memcpy(&dsa->ja[1], &ja[1], dsa->ne * sizeof(int));
         ufree(ja);
         dsa->ar = ucalloc(1+dsa->size, sizeof(double));
         memcpy(&dsa->ar[1], &ar[1], dsa->ne * sizeof(double));
         ufree(ar);
      }
      insist(dsa->ne < dsa->size);
      dsa->ne++;
      dsa->ia[dsa->ne] = i;
      dsa->ja[dsa->ne] = j;
      dsa->ar[dsa->ne] = aij;
      return;
}

static void transform(struct dsa *dsa)
{     /* transform original LP to standard formulation */
      LPX *lp = dsa->lp;
      int orig_m = dsa->orig_m;
      int orig_n = dsa->orig_n;
      int *ref = dsa->ref;
      int m = dsa->m;
      int n = dsa->n;
      double *b = dsa->b;
      double *c = dsa->c;
      int i, j, k, type, t, ii, len, *ind;
      double lb, ub, coef, rii, sjj, *val;
      /* initialize components of transformed LP */
      dsa->ne = 0;
      for (i = 1; i <= m; i++) b[i] = 0.0;
      c[0] = lpx_get_obj_coef(lp, 0);
      for (j = 1; j <= n; j++) c[j] = 0.0;
      /* i and j are, respectively, ordinal number of current row and
         ordinal number of current column in transformed LP */
      i = j = 0;
      /* transform rows (auxiliary variables) */
      for (k = 1; k <= orig_m; k++)
      {  type = lpx_get_row_type(lp, k);
         rii = lpx_get_rii(lp, k);
         lb = lpx_get_row_lb(lp, k) * rii;
         ub = lpx_get_row_ub(lp, k) * rii;
         switch (type)
         {  case LPX_FR:
               /* source: -inf < (L.F.) < +inf */
               /* result: ignore free row */
               ref[k] = 0;
               break;
            case LPX_LO:
               /* source: lb <= (L.F.) < +inf */
               /* result: (L.F.) - x' = lb, x' >= 0 */
               i++; j++;
               ref[k] = i;
               new_coef(dsa, i, j, -1.0);
               b[i] = lb;
               break;
            case LPX_UP:
               /* source: -inf < (L.F.) <= ub */
               /* result: (L.F.) + x' = ub, x' >= 0 */
               i++; j++;
               ref[k] = i;
               new_coef(dsa, i, j, +1.0);
               b[i] = ub;
               break;
            case LPX_DB:
               /* source: lb <= (L.F.) <= ub */
               /* result: (L.F.) - x' = lb, x' + x'' = ub - lb */
               i++; j++;
               ref[k] = i;
               new_coef(dsa, i, j, -1.0);
               b[i] = lb;
               i++;
               new_coef(dsa, i, j, +1.0);
               j++;
               new_coef(dsa, i, j, +1.0);
               b[i] = ub - lb;
               break;
            case LPX_FX:
               /* source: (L.F.) = lb */
               /* result: (L.F.) = lb */
               i++;
               ref[k] = i;
               b[i] = lb;
               break;
            default:
               insist(type != type);
         }
      }
      /* transform columns (structural variables) */
      ind = ucalloc(1+orig_m, sizeof(int));
      val = ucalloc(1+orig_m, sizeof(double));
      for (k = 1; k <= orig_n; k++)
      {  type = lpx_get_col_type(lp, k);
         sjj = lpx_get_sjj(lp, k);
         lb = lpx_get_col_lb(lp, k) / sjj;
         ub = lpx_get_col_ub(lp, k) / sjj;
         coef = lpx_get_obj_coef(lp, k) * sjj;
         len = lpx_get_mat_col(lp, k, ind, val);
         for (t = 1; t <= len; t++)
            val[t] *= (lpx_get_rii(lp, ind[t]) * sjj);
         switch (type)
         {  case LPX_FR:
               /* source: -inf < x < +inf */
               /* result: x = x' - x'', x' >= 0, x'' >= 0 */
               j++;
               ref[orig_m+k] = j;
               for (t = 1; t <= len; t++)
               {  ii = ref[ind[t]];
                  if (ii != 0) new_coef(dsa, ii, j, +val[t]);
               }
               c[j] = +coef;
               j++;
               for (t = 1; t <= len; t++)
               {  ii = ref[ind[t]];
                  if (ii != 0) new_coef(dsa, ii, j, -val[t]);
               }
               c[j] = -coef;
               break;
            case LPX_LO:
               /* source: lb <= x < +inf */
               /* result: x = lb + x', x' >= 0 */
               j++;
               ref[orig_m+k] = j;
               for (t = 1; t <= len; t++)
               {  ii = ref[ind[t]];
                  if (ii != 0)
                  {  new_coef(dsa, ii, j, val[t]);
                     b[ii] -= val[t] * lb;
                  }
               }
               c[j] = +coef;
               c[0] += c[j] * lb;
               break;
            case LPX_UP:
               /* source: -inf < x <= ub */
               /* result: x = ub - x', x' >= 0 */
               j++;
               ref[orig_m+k] = j;
               for (t = 1; t <= len; t++)
               {  ii = ref[ind[t]];
                  if (ii != 0)
                  {  new_coef(dsa, ii, j, -val[t]);
                     b[ii] -= val[t] * ub;
                  }
               }
               c[j] = -coef;
               c[0] -= c[j] * ub;
               break;
            case LPX_DB:
               /* source: lb <= x <= ub */
               /* result: x = lb + x', x' + x'' = ub - lb */
               j++;
               ref[orig_m+k] = j;
               for (t = 1; t <= len; t++)
               {  ii = ref[ind[t]];
                  if (ii != 0)
                  {  new_coef(dsa, ii, j, val[t]);
                     b[ii] -= val[t] * lb;
                  }
               }
               c[j] = +coef;
               c[0] += c[j] * lb;
               i++;
               new_coef(dsa, i, j, +1.0);
               j++;
               new_coef(dsa, i, j, +1.0);
               b[i] = ub - lb;
               break;
            case LPX_FX:
               /* source: x = lb */
               /* result: just substitute */
               ref[orig_m+k] = 0;
               for (t = 1; t <= len; t++)
               {  ii = ref[ind[t]];
                  if (ii != 0) b[ii] -= val[t] * lb;
               }
               c[0] += coef * lb;
               break;
            default:
               insist(type != type);
         }
      }
      ufree(ind);
      ufree(val);
      /* end of transformation */
      insist(i == m && j == n);
      /* change the objective sign in case of maximization */
      if (lpx_get_obj_dir(lp) == LPX_MAX)
         for (j = 0; j <= n; j++) c[j] = -c[j];
      return;
}

static void restore(struct dsa *dsa, double row_pval[],
      double row_dval[], double col_pval[], double col_dval[])
{     /* restore solution of original LP */
      LPX *lp = dsa->lp;
      int orig_m = dsa->orig_m;
      int orig_n = dsa->orig_n;
      int *ref = dsa->ref;
      int m = dsa->m;
      double *x = dsa->x;
      double *y = dsa->y;
      int dir = lpx_get_obj_dir(lp);
      int i, j, k, type, t, len, *ind;
      double lb, ub, rii, sjj, temp, *val;
      /* compute primal values of structural variables */
      for (k = 1; k <= orig_n; k++)
      {  j = ref[orig_m+k];
         type = lpx_get_col_type(lp, k);
         sjj = lpx_get_sjj(lp, k);
         lb = lpx_get_col_lb(lp, k) / sjj;
         ub = lpx_get_col_ub(lp, k) / sjj;
         switch (type)
         {  case LPX_FR:
               /* source: -inf < x < +inf */
               /* result: x = x' - x'', x' >= 0, x'' >= 0 */
               col_pval[k] = x[j] - x[j+1];
               break;
            case LPX_LO:
               /* source: lb <= x < +inf */
               /* result: x = lb + x', x' >= 0 */
               col_pval[k] = lb + x[j];
               break;
            case LPX_UP:
               /* source: -inf < x <= ub */
               /* result: x = ub - x', x' >= 0 */
               col_pval[k] = ub - x[j];
               break;
            case LPX_DB:
               /* source: lb <= x <= ub */
               /* result: x = lb + x', x' + x'' = ub - lb */
               col_pval[k] = lb + x[j];
               break;
            case LPX_FX:
               /* source: x = lb */
               /* result: just substitute */
               col_pval[k] = lb;
               break;
            default:
               insist(type != type);
         }
      }
      /* compute primal values of auxiliary variables */
      /* xR = A * xS */
      ind = ucalloc(1+orig_n, sizeof(int));
      val = ucalloc(1+orig_n, sizeof(double));
      for (k = 1; k <= orig_m; k++)
      {  rii = lpx_get_rii(lp, k);
         temp = 0.0;
         len = lpx_get_mat_row(lp, k, ind, val);
         for (t = 1; t <= len; t++)
         {  sjj = lpx_get_sjj(lp, ind[t]);
            temp += (rii * val[t] * sjj) * col_pval[ind[t]];
         }
         row_pval[k] = temp;
      }
      ufree(ind);
      ufree(val);
      /* compute dual values of auxiliary variables */
      for (k = 1; k <= orig_m; k++)
      {  type = lpx_get_row_type(lp, k);
         i = ref[k];
         switch (type)
         {  case LPX_FR:
               insist(i == 0);
               row_dval[k] = 0.0;
               break;
            case LPX_LO:
            case LPX_UP:
            case LPX_DB:
            case LPX_FX:
               insist(1 <= i && i <= m);
               row_dval[k] = (dir == LPX_MIN ? +1.0 : -1.0) * y[i];
               break;
            default:
               insist(type != type);
         }
      }
      /* compute dual values of structural variables */
      /* dS = cS - A' * (dR - cR) */
      ind = ucalloc(1+orig_m, sizeof(int));
      val = ucalloc(1+orig_m, sizeof(double));
      for (k = 1; k <= orig_n; k++)
      {  sjj = lpx_get_sjj(lp, k);
         temp = lpx_get_obj_coef(lp, k) / sjj;
         len = lpx_get_mat_col(lp, k, ind, val);
         for (t = 1; t <= len; t++)
         {  rii = lpx_get_rii(lp, ind[t]);
            temp -= (rii * val[t] * sjj) * row_dval[ind[t]];
         }
         col_dval[k] = temp;
      }
      ufree(ind);
      ufree(val);
      /* unscale solution of original LP */
      for (i = 1; i <= orig_m; i++)
      {  rii = lpx_get_rii(lp, i);
         row_pval[i] /= rii;
         row_dval[i] *= rii;
      }
      for (j = 1; j <= orig_n; j++)
      {  sjj = lpx_get_sjj(lp, j);
         col_pval[j] *= sjj;
         col_dval[j] /= sjj;
      }
      return;
}

int lpx_interior(LPX *_lp)
{     /* easy-to-use driver to the interior-point method */
      struct dsa _dsa, *dsa = &_dsa;
      int ret, status;
      double *row_pval, *row_dval, *col_pval, *col_dval;
      /* initialize working area */
      dsa->lp = _lp;
      dsa->orig_m = lpx_get_num_rows(dsa->lp);
      dsa->orig_n = lpx_get_num_cols(dsa->lp);
      dsa->ref = NULL;
      dsa->m = 0;
      dsa->n = 0;
      dsa->size = 0;
      dsa->ne = 0;
      dsa->ia = NULL;
      dsa->ja = NULL;
      dsa->ar = NULL;
      dsa->b = NULL;
      dsa->c = NULL;
      dsa->x = NULL;
      dsa->y = NULL;
      dsa->z = NULL;
      /* check if the problem is empty */
      if (!(dsa->orig_m > 0 && dsa->orig_n > 0))
      {  print("lpx_interior: problem has no rows and/or columns");
         ret = LPX_E_FAULT;
         goto done;
      }
      /* issue warning about dense columns */
      if (dsa->orig_m >= 200)
      {  int j, len, ndc = 0;
         for (j = 1; j <= dsa->orig_n; j++)
         {  len = lpx_get_mat_col(dsa->lp, j, NULL, NULL);
            if ((double)len > 0.30 * (double)dsa->orig_m) ndc++;
         }
         if (ndc == 1)
            print("lpx_interior: WARNING: PROBLEM HAS ONE DENSE COLUMN")
               ;
         else if (ndc > 0)
            print("lpx_interior: WARNING: PROBLEM HAS %d DENSE COLUMNS",
               ndc);
      }
      /* determine dimension of transformed LP */
      print("lpx_interior: original LP problem has %d rows and %d colum"
         "ns", dsa->orig_m, dsa->orig_n);
      calc_mn(dsa);
      print("lpx_interior: transformed LP problem has %d rows and %d co"
         "lumns", dsa->m, dsa->n);
      /* transform original LP to standard formulation */
      dsa->ref = ucalloc(1+dsa->orig_m+dsa->orig_n, sizeof(int));
      dsa->size = lpx_get_num_nz(dsa->lp);
      if (dsa->size < dsa->m) dsa->size = dsa->m;
      if (dsa->size < dsa->n) dsa->size = dsa->n;
      dsa->ne = 0;
      dsa->ia = ucalloc(1+dsa->size, sizeof(int));
      dsa->ja = ucalloc(1+dsa->size, sizeof(int));
      dsa->ar = ucalloc(1+dsa->size, sizeof(double));
      dsa->b = ucalloc(1+dsa->m, sizeof(double));
      dsa->c = ucalloc(1+dsa->n, sizeof(double));
      transform(dsa);
      /* solve the transformed LP problem */
      {  int *A_ptr = ucalloc(1+dsa->m+1, sizeof(int));
         int *A_ind = ucalloc(1+dsa->ne, sizeof(int));
         double *A_val = ucalloc(1+dsa->ne, sizeof(double));
         int i, k, pos, len;
         /* determine row lengths */
         for (i = 1; i <= dsa->m; i++)
            A_ptr[i] = 0;
         for (k = 1; k <= dsa->ne; k++)
            A_ptr[dsa->ia[k]]++;
         /* set up row pointers */
         pos = 1;
         for (i = 1; i <= dsa->m; i++)
            len = A_ptr[i], pos += len, A_ptr[i] = pos;
         A_ptr[dsa->m+1] = pos;
         insist(pos - 1 == dsa->ne);
         /* build the matrix */
         for (k = 1; k <= dsa->ne; k++)
         {  pos = --A_ptr[dsa->ia[k]];
            A_ind[pos] = dsa->ja[k];
            A_val[pos] = dsa->ar[k];
         }
         ufree(dsa->ia), dsa->ia = NULL;
         ufree(dsa->ja), dsa->ja = NULL;
         ufree(dsa->ar), dsa->ar = NULL;
         dsa->x = ucalloc(1+dsa->n, sizeof(double));
         dsa->y = ucalloc(1+dsa->m, sizeof(double));
         dsa->z = ucalloc(1+dsa->n, sizeof(double));
         ret = ipm_main(dsa->m, dsa->n, A_ptr, A_ind, A_val, dsa->b,
            dsa->c, dsa->x, dsa->y, dsa->z);
         ufree(A_ptr);
         ufree(A_ind);
         ufree(A_val);
         ufree(dsa->b), dsa->b = NULL;
         ufree(dsa->c), dsa->c = NULL;
      }
      /* analyze return code reported by the solver */
      switch (ret)
      {  case 0:
            /* optimal solution found */
            ret = LPX_E_OK;
            break;
         case 1:
            /* problem has no feasible (primal or dual) solution */
            ret = LPX_E_NOFEAS;
            break;
         case 2:
            /* no convergence */
            ret = LPX_E_NOCONV;
            break;
         case 3:
            /* iterations limit exceeded */
            ret = LPX_E_ITLIM;
            break;
         case 4:
            /* numerical instability on solving Newtonian system */
            ret = LPX_E_INSTAB;
            break;
         default:
            insist(ret != ret);
      }
      /* recover solution of original LP */
      row_pval = ucalloc(1+dsa->orig_m, sizeof(double));
      row_dval = ucalloc(1+dsa->orig_m, sizeof(double));
      col_pval = ucalloc(1+dsa->orig_n, sizeof(double));
      col_dval = ucalloc(1+dsa->orig_n, sizeof(double));
      restore(dsa, row_pval, row_dval, col_pval, col_dval);
      ufree(dsa->ref), dsa->ref = NULL;
      ufree(dsa->x), dsa->x = NULL;
      ufree(dsa->y), dsa->y = NULL;
      ufree(dsa->z), dsa->z = NULL;
      /* store solution components into the problem object */
      status = (ret == LPX_E_OK ? LPX_T_OPT : LPX_T_UNDEF);
      lpx_put_ipt_soln(dsa->lp, status, row_pval, row_dval, col_pval,
         col_dval);
      ufree(row_pval);
      ufree(row_dval);
      ufree(col_pval);
      ufree(col_dval);
done: /* return to the calling program */
      return ret;
}

/* eof */
