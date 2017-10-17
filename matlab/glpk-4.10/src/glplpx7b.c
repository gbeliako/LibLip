/* glplpx7b.c */

/* (reserved for copyright notice) */

#include <float.h>
#include <limits.h>
#include <string.h>
#include "glplib.h"
#include "glplpx.h"

static double get_row_lb(LPX *lp, int i)
{     /* this routine returns lower bound of row i or -DBL_MAX if the
         row has no lower bound */
      double lb;
      switch (lpx_get_row_type(lp, i))
      {  case LPX_FR:
         case LPX_UP:
            lb = -DBL_MAX;
            break;
         case LPX_LO:
         case LPX_DB:
         case LPX_FX:
            lb = lpx_get_row_lb(lp, i);
            break;
         default:
            insist(lp != lp);
      }
      return lb;
}

static double get_row_ub(LPX *lp, int i)
{     /* this routine returns upper bound of row i or +DBL_MAX if the
         row has no upper bound */
      double ub;
      switch (lpx_get_row_type(lp, i))
      {  case LPX_FR:
         case LPX_LO:
            ub = +DBL_MAX;
            break;
         case LPX_UP:
         case LPX_DB:
         case LPX_FX:
            ub = lpx_get_row_ub(lp, i);
            break;
         default:
            insist(lp != lp);
      }
      return ub;
}

static double get_col_lb(LPX *lp, int j)
{     /* this routine returns lower bound of column j or -DBL_MAX if
         the column has no lower bound */
      double lb;
      switch (lpx_get_col_type(lp, j))
      {  case LPX_FR:
         case LPX_UP:
            lb = -DBL_MAX;
            break;
         case LPX_LO:
         case LPX_DB:
         case LPX_FX:
            lb = lpx_get_col_lb(lp, j);
            break;
         default:
            insist(lp != lp);
      }
      return lb;
}

static double get_col_ub(LPX *lp, int j)
{     /* this routine returns upper bound of column j or +DBL_MAX if
         the column has no upper bound */
      double ub;
      switch (lpx_get_col_type(lp, j))
      {  case LPX_FR:
         case LPX_LO:
            ub = +DBL_MAX;
            break;
         case LPX_UP:
         case LPX_DB:
         case LPX_FX:
            ub = lpx_get_col_ub(lp, j);
            break;
         default:
            insist(lp != lp);
      }
      return ub;
}

static int is_binary(LPX *lp, int j)
{     /* this routine checks if variable x[j] is binary */
      return
         lpx_get_col_kind(lp, j) == LPX_IV &&
         lpx_get_col_type(lp, j) == LPX_DB &&
         lpx_get_col_lb(lp, j) == 0.0 && lpx_get_col_ub(lp, j) == 1.0;
}

static double eval_lf_min(LPX *lp, int len, int ind[], double val[])
{     /* this routine computes the minimum of a specified linear form

            sum a[j]*x[j]
             j

         using the formula:

            min =   sum   a[j]*lb[j] +   sum   a[j]*ub[j],
                  j in J+              j in J-

         where J+ = {j: a[j] > 0}, J- = {j: a[j] < 0}, lb[j] and ub[j]
         are lower and upper bound of variable x[j], resp. */
      int j, t;
      double lb, ub, sum;
      sum = 0.0;
      for (t = 1; t <= len; t++)
      {  j = ind[t];
         if (val[t] > 0.0)
         {  lb = get_col_lb(lp, j);
            if (lb == -DBL_MAX)
            {  sum = -DBL_MAX;
               break;
            }
            sum += val[t] * lb;
         }
         else if (val[t] < 0.0)
         {  ub = get_col_ub(lp, j);
            if (ub == +DBL_MAX)
            {  sum = -DBL_MAX;
               break;
            }
            sum += val[t] * ub;
         }
         else
            insist(val != val);
      }
      return sum;
}

static double eval_lf_max(LPX *lp, int len, int ind[], double val[])
{     /* this routine computes the maximum of a specified linear form

            sum a[j]*x[j]
             j

         using the formula:

            max =   sum   a[j]*ub[j] +   sum   a[j]*lb[j],
                  j in J+              j in J-

         where J+ = {j: a[j] > 0}, J- = {j: a[j] < 0}, lb[j] and ub[j]
         are lower and upper bound of variable x[j], resp. */
      int j, t;
      double lb, ub, sum;
      sum = 0.0;
      for (t = 1; t <= len; t++)
      {  j = ind[t];
         if (val[t] > 0.0)
         {  ub = get_col_ub(lp, j);
            if (ub == +DBL_MAX)
            {  sum = +DBL_MAX;
               break;
            }
            sum += val[t] * ub;
         }
         else if (val[t] < 0.0)
         {  lb = get_col_lb(lp, j);
            if (lb == -DBL_MAX)
            {  sum = +DBL_MAX;
               break;
            }
            sum += val[t] * lb;
         }
         else
            insist(val != val);
      }
      return sum;
}

/*----------------------------------------------------------------------
-- probing - determine logical relation between binary variables.
--
-- This routine tentatively sets a binary variable to 0 and then to 1
-- and examines whether another binary variable is caused to be fixed.
--
-- The examination is based only on one row (constraint), which is the
-- following:
--
--    L <= sum a[j]*x[j] <= U.                                       (1)
--          j
--
-- Let x[p] be a probing variable, x[q] be an examined variable. Then
-- (1) can be written as:
--
--    L <=   sum  a[j]*x[j] + a[p]*x[p] + a[q]*x[q] <= U,            (2)
--         j in J'
--
-- where J' = {j: j != p and j != q}.
--
-- Let
--
--    L' = L - a[p]*x[p],                                            (3)
--
--    U' = U - a[p]*x[p],                                            (4)
--
-- where x[p] is assumed to be fixed at 0 or 1. So (2) can be rewritten
-- as follows:
--
--    L' <=   sum  a[j]*x[j] + a[q]*x[q] <= U',                      (5)
--          j in J'
--
-- from where we have:
--
--    L' -  sum  a[j]*x[j] <= a[q]*x[q] <= U' -  sum  a[j]*x[j].     (6)
--        j in J'                              j in J'
--
-- Thus,
--
--    min a[q]*x[q] = L' - MAX,                                      (7)
--
--    max a[q]*x[q] = U' - MIN,                                      (8)
--
-- where
--
--    MIN = min  sum  a[j]*x[j],                                     (9)
--             j in J'
--
--    MAX = max  sum  a[j]*x[j].                                    (10)
--             j in J'
--
-- Formulae (7) and (8) allows determining implied lower and upper
-- bounds of x[q].
--
-- Parameters len, val, L and U specify the constraint (1).
--
-- Parameters lf_min and lf_max specify implied lower and upper bounds
-- of the linear form (1). It is assumed that these bounds are computed
-- with the routines eval_lf_min and eval_lf_max (see above).
--
-- Parameter p specifies the probing variable x[p], which is set to 0
-- (if set is 0) or to 1 (if set is 1).
--
-- Parameter q specifies the examined variable x[q].
--
-- On exit the routine returns one of the following codes:
--
-- 0 - there is no logical relation between x[p] and x[q];
-- 1 - x[q] can take only on value 0;
-- 2 - x[q] can take only on value 1. */

static int probing(int len, double val[], double L, double U,
      double lf_min, double lf_max, int p, int set, int q)
{     double temp;
      insist(1 <= p && p < q && q <= len);
      /* compute L' (3) */
      if (L != -DBL_MAX && set) L -= val[p];
      /* compute U' (4) */
      if (U != +DBL_MAX && set) U -= val[p];
      /* compute MIN (9) */
      if (lf_min != -DBL_MAX)
      {  if (val[p] < 0.0) lf_min -= val[p];
         if (val[q] < 0.0) lf_min -= val[q];
      }
      /* compute MAX (10) */
      if (lf_max != +DBL_MAX)
      {  if (val[p] > 0.0) lf_max -= val[p];
         if (val[q] > 0.0) lf_max -= val[q];
      }
      /* compute implied lower bound of x[q]; see (7), (8) */
      if (val[q] > 0.0)
      {  if (L == -DBL_MAX || lf_max == +DBL_MAX)
            temp = -DBL_MAX;
         else
            temp = (L - lf_max) / val[q];
      }
      else
      {  if (U == +DBL_MAX || lf_min == -DBL_MAX)
            temp = -DBL_MAX;
         else
            temp = (U - lf_min) / val[q];
      }
      if (temp > 0.001) return 2;
      /* compute implied upper bound of x[q]; see (7), (8) */
      if (val[q] > 0.0)
      {  if (U == +DBL_MAX || lf_min == -DBL_MAX)
            temp = +DBL_MAX;
         else
            temp = (U - lf_min) / val[q];
      }
      else
      {  if (L == -DBL_MAX || lf_max == +DBL_MAX)
            temp = +DBL_MAX;
         else
            temp = (L - lf_max) / val[q];
      }
      if (temp < 0.999) return 1;
      /* there is no logical relation between x[p] and x[q] */
      return 0;
}

struct COG
{     /* conflict graph; it represents logical relations between binary
         variables and has a vertex for each binary variable and its
         complement, and an edge between two vertices when at most one
         of the variables represented by the vertices can equal one in
         an optimal solution */
      int n;
      /* number of variables */
      int nb;
      /* number of binary variables represented in the graph (note that
         not all binary variables can be represented); vertices which
         correspond to binary variables have numbers 1, ..., nb while
         vertices which correspond to complements of binary variables
         have numbers nb+1, ..., nb+nb */
      int ne;
      /* number of edges in the graph */
      int *vert; /* int vert[1+n]; */
      /* if x[j] is a binary variable represented in the graph, vert[j]
         is the vertex number corresponding to x[j]; otherwise vert[j]
         is zero */
      int *orig; /* int list[1:nb]; */
      /* if vert[j] = k > 0, then orig[k] = j */
      unsigned char *a;
      /* adjacency matrix of the graph having 2*nb rows and columns;
         only strict lower triangle is stored in dense packed form */
};

/*----------------------------------------------------------------------
-- lpx_create_cog - create the conflict graph.
--
-- SYNOPSIS
--
-- #include "glplpx.h"
-- void *lpx_create_cog(LPX *lp);
--
-- DESCRIPTION
--
-- The routine lpx_create_cog creates the conflict graph for a given
-- problem instance.
--
-- RETURNS
--
-- If the graph has been created, the routine returns a pointer to it.
-- Otherwise the routine returns NULL. */

#define MAX_NB 4000
#define MAX_ROW_LEN 500

void *lpx_create_cog(LPX *lp)
{     struct COG *cog = NULL;
      int m, n, nb, i, j, p, q, len, *ind, *vert, *orig;
      double L, U, lf_min, lf_max, *val;
      print("Creating the conflict graph...");
      m = lpx_get_num_rows(lp);
      n = lpx_get_num_cols(lp);
      /* determine which binary variables should be included in the
         conflict graph */
      nb = 0;
      vert = ucalloc(1+n, sizeof(int));
      for (j = 1; j <= n; j++) vert[j] = 0;
      orig = ucalloc(1+n, sizeof(int));
      ind = ucalloc(1+n, sizeof(int));
      val = ucalloc(1+n, sizeof(double));
      for (i = 1; i <= m; i++)
      {  L = get_row_lb(lp, i);
         U = get_row_ub(lp, i);
         if (L == -DBL_MAX && U == +DBL_MAX) continue;
         len = lpx_get_mat_row(lp, i, ind, val);
         if (len > MAX_ROW_LEN) continue;
         lf_min = eval_lf_min(lp, len, ind, val);
         lf_max = eval_lf_max(lp, len, ind, val);
         for (p = 1; p <= len; p++)
         {  if (!is_binary(lp, ind[p])) continue;
            for (q = p+1; q <= len; q++)
            {  if (!is_binary(lp, ind[q])) continue;
               if (probing(len, val, L, U, lf_min, lf_max, p, 0, q) ||
                   probing(len, val, L, U, lf_min, lf_max, p, 1, q))
               {  /* there is a logical relation */
                  /* include the first variable in the graph */
                  j = ind[p];
                  if (vert[j] == 0) nb++, vert[j] = nb, orig[nb] = j;
                  /* incude the second variable in the graph */
                  j = ind[q];
                  if (vert[j] == 0) nb++, vert[j] = nb, orig[nb] = j;
               }
            }
         }
      }
      /* if the graph is either empty or has too many vertices, do not
         create it */
      if (nb == 0 || nb > MAX_NB)
      {  print("The conflict graph is either empty or too big");
         ufree(vert);
         ufree(orig);
         goto done;
      }
      /* create the conflict graph */
      cog = umalloc(sizeof(struct COG));
      cog->n = n;
      cog->nb = nb;
      cog->ne = 0;
      cog->vert = vert;
      cog->orig = orig;
      len = nb + nb; /* number of vertices */
      len = (len * (len - 1)) / 2; /* number of entries in triangle */
      len = (len + (CHAR_BIT - 1)) / CHAR_BIT; /* bytes needed */
      cog->a = umalloc(len);
      memset(cog->a, 0, len);
      for (j = 1; j <= nb; j++)
      {  /* add edge between variable and its complement */
         lpx_add_cog_edge(cog, +orig[j], -orig[j]);
      }
      for (i = 1; i <= m; i++)
      {  L = get_row_lb(lp, i);
         U = get_row_ub(lp, i);
         if (L == -DBL_MAX && U == +DBL_MAX) continue;
         len = lpx_get_mat_row(lp, i, ind, val);
         if (len > MAX_ROW_LEN) continue;
         lf_min = eval_lf_min(lp, len, ind, val);
         lf_max = eval_lf_max(lp, len, ind, val);
         for (p = 1; p <= len; p++)
         {  if (!is_binary(lp, ind[p])) continue;
            for (q = p+1; q <= len; q++)
            {  if (!is_binary(lp, ind[q])) continue;
               /* set x[p] to 0 and examine x[q] */
               switch (probing(len, val, L, U, lf_min, lf_max, p, 0, q))
               {  case 0:
                     /* no logical relation */
                     break;
                  case 1:
                     /* x[p] = 0 implies x[q] = 0 */
                     lpx_add_cog_edge(cog, -ind[p], +ind[q]);
                     break;
                  case 2:
                     /* x[p] = 0 implies x[q] = 1 */
                     lpx_add_cog_edge(cog, -ind[p], -ind[q]);
                     break;
                  default:
                     insist(lp != lp);
               }
               /* set x[p] to 1 and examine x[q] */
               switch (probing(len, val, L, U, lf_min, lf_max, p, 1, q))
               {  case 0:
                     /* no logical relation */
                     break;
                  case 1:
                     /* x[p] = 1 implies x[q] = 0 */
                     lpx_add_cog_edge(cog, +ind[p], +ind[q]);
                     break;
                  case 2:
                     /* x[p] = 1 implies x[q] = 1 */
                     lpx_add_cog_edge(cog, +ind[p], -ind[q]);
                     break;
                  default:
                     insist(lp != lp);
               }
            }
         }
      }
      print("The conflict graph has 2*%d vertices and %d edges",
         cog->nb, cog->ne);
done: ufree(ind);
      ufree(val);
      return cog;
}

/*----------------------------------------------------------------------
-- lpx_add_cog_edge - add edge to the conflict graph.
--
-- SYNOPSIS
--
-- #include "glplpx.h"
-- void lpx_add_cog_edge(void *cog, int i, int j);
--
-- DESCRIPTION
--
-- The routine lpx_add_cog_edge adds an edge to the conflict graph.
-- The edge connects x[i] (if i > 0) or its complement (if i < 0) and
-- x[j] (if j > 0) or its complement (if j < 0), where i and j are
-- original ordinal numbers of corresponding variables. */

void lpx_add_cog_edge(void *_cog, int i, int j)
{     struct COG *cog = _cog;
      int k;
      insist(i != j);
      /* determine indices of corresponding vertices */
      if (i > 0)
      {  insist(1 <= i && i <= cog->n);
         i = cog->vert[i];
         insist(i != 0);
      }
      else
      {  i = -i;
         insist(1 <= i && i <= cog->n);
         i = cog->vert[i];
         insist(i != 0);
         i += cog->nb;
      }
      if (j > 0)
      {  insist(1 <= j && j <= cog->n);
         j = cog->vert[j];
         insist(j != 0);
      }
      else
      {  j = -j;
         insist(1 <= j && j <= cog->n);
         j = cog->vert[j];
         insist(j != 0);
         j += cog->nb;
      }
      /* only lower triangle is stored, so we need i > j */
      if (i < j) k = i, i = j, j = k;
      k = ((i - 1) * (i - 2)) / 2 + (j - 1);
      cog->a[k / CHAR_BIT] |=
         (unsigned char)(1 << ((CHAR_BIT - 1) - k % CHAR_BIT));
      cog->ne++;
      return;
}

/*----------------------------------------------------------------------
-- MAXIMUM WEIGHT CLIQUE
--
-- Two subroutines sub() and wclique() below are intended to find a
-- maximum weight clique in a given undirected graph. These subroutines
-- are slightly modified version of the program WCLIQUE developed by
-- Patric Ostergard <http://www.tcs.hut.fi/~pat/wclique.html> and based
-- on ideas from the article "P. R. J. Ostergard, A new algorithm for
-- the maximum-weight clique problem, submitted for publication", which
-- in turn is a generalization of the algorithm for unweighted graphs
-- presented in "P. R. J. Ostergard, A fast algorithm for the maximum
-- clique problem, submitted for publication".
--
-- USED WITH PERMISSION OF THE AUTHOR OF THE ORIGINAL CODE. */

struct dsa
{     /* dynamic storage area */
      int n;
      /* number of vertices */
      int *wt; /* int wt[0:n-1]; */
      /* weights */
      unsigned char *a;
      /* adjacency matrix (packed lower triangle without main diag.) */
      int record;
      /* weight of best clique */
      int rec_level;
      /* number of vertices in best clique */
      int *rec; /* int rec[0:n-1]; */
      /* best clique so far */
      int *clique; /* int clique[0:n-1]; */
      /* table for pruning */
      int *set; /* int set[0:n-1]; */
      /* current clique */
};

#define n         (dsa->n)
#define wt        (dsa->wt)
#define a         (dsa->a)
#define record    (dsa->record)
#define rec_level (dsa->rec_level)
#define rec       (dsa->rec)
#define clique    (dsa->clique)
#define set       (dsa->set)

#if 0
static int is_edge(struct dsa *dsa, int i, int j)
{     /* if there is arc (i,j), the routine returns true; otherwise
         false; 0 <= i, j < n */
      int k;
      insist(0 <= i && i < n);
      insist(0 <= j && j < n);
      if (i == j) return 0;
      if (i < j) k = i, i = j, j = k;
      k = (i * (i - 1)) / 2 + j;
      return a[k / CHAR_BIT] &
         (unsigned char)(1 << ((CHAR_BIT - 1) - k % CHAR_BIT));
}
#else
#define is_edge(dsa, i, j) ((i) == (j) ? 0 : \
      (i) > (j) ? is_edge1(i, j) : is_edge1(j, i))
#define is_edge1(i, j) is_edge2(((i) * ((i) - 1)) / 2 + (j))
#define is_edge2(k) (a[(k) / CHAR_BIT] & \
      (unsigned char)(1 << ((CHAR_BIT - 1) - (k) % CHAR_BIT)))
#endif

static void sub(struct dsa *dsa, int ct, int table[], int level,
      int weight, int l_weight)
{     int i, j, k, curr_weight, left_weight, *p1, *p2, *newtable;
      newtable = ucalloc(n, sizeof(int));
      if (ct <= 0)
      {  /* 0 or 1 elements left; include these */
         if (ct == 0)
         {  set[level++] = table[0];
            weight += l_weight;
         }
         if (weight > record)
         {  record = weight;
            rec_level = level;
            for (i = 0; i < level; i++) rec[i] = set[i];
         }
         goto done;
      }
      for (i = ct; i >= 0; i--)
      {  if ((level == 0) && (i < ct)) goto done;
         k = table[i];
         if ((level > 0) && (clique[k] <= (record - weight)))
            goto done; /* prune */
         set[level] = k;
         curr_weight = weight + wt[k];
         l_weight -= wt[k];
         if (l_weight <= (record - curr_weight))
            goto done; /* prune */
         p1 = newtable;
         p2 = table;
         left_weight = 0;
         while (p2 < table + i)
         {  j = *p2++;
            if (is_edge(dsa, j, k))
            {  *p1++ = j;
               left_weight += wt[j];
            }
         }
         if (left_weight <= (record - curr_weight)) continue;
         sub(dsa, p1 - newtable - 1, newtable, level + 1, curr_weight,
            left_weight);
      }
done: ufree(newtable);
      return;
}

static int wclique(int _n, int w[], unsigned char _a[], int sol[])
{     struct dsa _dsa, *dsa = &_dsa;
      int i, j, p, max_wt, max_nwt, wth, *used, *nwt, *pos;
      double timer;
      n = _n;
      wt = &w[1];
      a = _a;
      record = 0;
      rec_level = 0;
      rec = &sol[1];
      clique = ucalloc(n, sizeof(int));
      set = ucalloc(n, sizeof(int));
      used = ucalloc(n, sizeof(int));
      nwt = ucalloc(n, sizeof(int));
      pos = ucalloc(n, sizeof(int));
      /* start timer */
      timer = utime();
      /* order vertices */
      for (i = 0; i < n; i++)
      {  nwt[i] = 0;
         for (j = 0; j < n; j++)
            if (is_edge(dsa, i, j)) nwt[i] += wt[j];
      }
      for (i = 0; i < n; i++)
         used[i] = 0;
      for (i = n-1; i >= 0; i--)
      {  max_wt = -1;
         max_nwt = -1;
         for (j = 0; j < n; j++)
         {  if ((!used[j]) && ((wt[j] > max_wt) || (wt[j] == max_wt
               && nwt[j] > max_nwt)))
            {  max_wt = wt[j];
               max_nwt = nwt[j];
               p = j;
            }
         }
         pos[i] = p;
         used[p] = 1;
         for (j = 0; j < n; j++)
            if ((!used[j]) && (j != p) && (is_edge(dsa, p, j)))
               nwt[j] -= wt[p];
      }
      /* main routine */
      wth = 0;
      for (i = 0; i < n; i++)
      {  wth += wt[pos[i]];
         sub(dsa, i, pos, 0, 0, wth);
         clique[pos[i]] = record;
         if (utime() >= timer + 5.0)
         {  /* print current record and reset timer */
            print("level = %d (%d); best = %d", i+1, n, record);
            timer = utime();
         }
      }
      ufree(clique);
      ufree(set);
      ufree(used);
      ufree(nwt);
      ufree(pos);
      /* return the solution found */
      for (i = 1; i <= rec_level; i++) sol[i]++;
      return rec_level;
}

#undef n
#undef wt
#undef a
#undef record
#undef rec_level
#undef rec
#undef clique
#undef set

/*----------------------------------------------------------------------
-- lpx_clique_cut - generate cluque cut.
--
-- SYNOPSIS
--
-- #include "glplpx.h"
-- int lpx_clique_cut(LPX *lp, void *cog, int ind[], double val[]);
--
-- DESCRIPTION
--
-- The routine lpx_clique_cut generates a clique cut using the conflict
-- graph specified by the parameter cog.
--
-- If a violated clique cut has been found, it has the following form:
--
--    sum{j in J} a[j]*x[j] <= b.
--
-- Variable indices j in J are stored in elements ind[1], ..., ind[len]
-- while corresponding constraint coefficients are stored in elements
-- val[1], ..., val[len], where len is returned on exit. The right-hand
-- side b is stored in element val[0].
--
-- RETURNS
--
-- If the cutting plane has been successfully generated, the routine
-- returns 1 <= len <= n, which is the number of non-zero coefficients
-- in the inequality constraint. Otherwise, the routine returns zero. */

int lpx_clique_cut(LPX *lp, void *_cog, int ind[], double val[])
{     struct COG *cog = _cog;
      int n = lpx_get_num_cols(lp);
      int j, t, v, card, temp, len = 0, *w, *sol;
      double x, sum, b, *vec;
      /* allocate working arrays */
      w = ucalloc(1 + 2 * cog->nb, sizeof(int));
      sol = ucalloc(1 + 2 * cog->nb, sizeof(int));
      vec = ucalloc(1+n, sizeof(double));
      /* assign weights to vertices of the conflict graph */
      for (t = 1; t <= cog->nb; t++)
      {  j = cog->orig[t];
         x = lpx_get_col_prim(lp, j);
         temp = (int)(100.0 * x + 0.5);
         if (temp < 0) temp = 0;
         if (temp > 100) temp = 100;
         w[t] = temp;
         w[cog->nb + t] = 100 - temp;
      }
      /* find a clique of maximum weight */
      card = wclique(2 * cog->nb, w, cog->a, sol);
      /* compute the clique weight for unscaled values */
      sum = 0.0;
      for ( t = 1; t <= card; t++)
      {  v = sol[t];
         insist(1 <= v && v <= 2 * cog->nb);
         if (v <= cog->nb)
         {  /* vertex v corresponds to binary variable x[j] */
            j = cog->orig[v];
            x = lpx_get_col_prim(lp, j);
            sum += x;
         }
         else
         {  /* vertex v corresponds to the complement of x[j] */
            j = cog->orig[v - cog->nb];
            x = lpx_get_col_prim(lp, j);
            sum += 1.0 - x;
         }
      }
      /* if the sum of binary variables and their complements in the
         clique greater than 1, the clique cut is violated */
      if (sum >= 1.01)
      {  /* construct the inquality */
         for (j = 1; j <= n; j++) vec[j] = 0;
         b = 1.0;
         for (t = 1; t <= card; t++)
         {  v = sol[t];
            if (v <= cog->nb)
            {  /* vertex v corresponds to binary variable x[j] */
               j = cog->orig[v];
               insist(1 <= j && j <= n);
               vec[j] += 1.0;
            }
            else
            {  /* vertex v corresponds to the complement of x[j] */
               j = cog->orig[v - cog->nb];
               insist(1 <= j && j <= n);
               vec[j] -= 1.0;
               b -= 1.0;
            }
         }
         insist(len == 0);
         for (j = 1; j <= n; j++)
         {  if (vec[j] != 0.0)
            {  len++;
               ind[len] = j, val[len] = vec[j];
            }
         }
         ind[0] = 0, val[0] = b;
      }
      /* free working arrays */
      ufree(w);
      ufree(sol);
      ufree(vec);
      /* return to the calling program */
      return len;
}

/*----------------------------------------------------------------------
-- lpx_delete_cog - delete the conflict graph.
--
-- SYNOPSIS
--
-- #include "glplpx.h"
-- void lpx_delete_cog(void *cog);
--
-- DESCRIPTION
--
-- The routine lpx_delete_cog deletes the conflict graph, which the
-- parameter cog points to, freeing all the memory allocated to this
-- object. */

void lpx_delete_cog(void *_cog)
{     struct COG *cog = _cog;
      ufree(cog->vert);
      ufree(cog->orig);
      ufree(cog->a);
      ufree(cog);
}

/* eof */
