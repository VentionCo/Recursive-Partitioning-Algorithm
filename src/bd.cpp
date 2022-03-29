/* Copyright (C) 2007-2016 Rafael Durbano Lobato
 *
 * This file is part of Recursive Partitioning Algorithm.
 *
 * Recursive Partitioning Algorithm is free software: you can
 * redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * Recursive Partitioning Algorithm is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Recursive Partitioning Algorithm. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * Contact information:
 *
 * Rafael Durbano Lobato <lobato@ime.usp.br>
 * http://www.ime.usp.br/~lobato/
 */

/* Implementation of the algorithm "Recursive Five-block Heuristic"
 * described in [1] for solving the Manufacturer's Pallet Loading
 * Problem.
 *
 * [1] E. G. Birgin, R. D. Lobato and R. Morabito. An effective
 *     recursive partitioning approach for the packing of identical
 *     rectangles in a rectangle. Journal of the Operational Research
 *     Society (2010) 61, 306-320.
 *
 * Last update: July 8th, 2009
 * Author:      Rafael Durbano Lobato
 * Contact:     lobato@ime.usp.br
 */

/******************************************************************
 ******************************************************************/

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/resource.h>
#include <sys/time.h>
#include <sys/times.h>

#include "bd.h"
#include "sets.h"
#include "util.h"

#define INFINITY_ 2000000000

#define lowerBound(L, W, l, w) std::max ((L / l) * (W / w), (L / w) * (W / l));

int BD (int L, int W, int l, int w, int n);

int barnesBound (int L, int W, int l, int w);

/* Maximum level of recursion (maximum tree search depth). */
int N = INFINITY_;

/* Arrays of indices for indexing the matrices that store
 * information about problems (L,W), where (L,W) belongs to X' x Y'
 * and X' and Y' are the raster points sets associated to (L,l,w) and
 * (W,l,w), respectively. */
extern int *indexX, *indexY;

/* Lower and upper bounds of each subproblem. */
extern int **lowerBound, **upperBound;

extern Set normalSetX;

/* Indicate in which level of the recursion (or in which depth of the
 * tree search) the solution was found. Initially, solutionDepth[L][W]
 * = N, for all (L,W) subproblem. If the optimal solution was found
 * for a subproblem (L,W), than solutionDepth[L][W] = -1. */
int **solutionDepth;

/* Array that stores the normalized values of each integer between 0 and L.
 * normalize[x] = max {r in X' | r <= x} */
extern int *normalize;

/* Store the points that determine the divisions of the rectangles. */
extern CutPoint **cutPoints;

/* Indicate if the limit of the recursion was reached during the
 * resolution of a problem. */
int **reachedLimit;

/******************************************************************
 ******************************************************************/

inline int
localUpperBound (int iX, int iY)
{
  if (solutionDepth[iX][iY] == -1)
    {
      return lowerBound[iX][iY];
    }
  else
    {
      return upperBound[iX][iY];
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Store the points x1, x2, y1 and y2 that determine the cut in the
 * rectangle (L,W).
 */
void
storeCutPoint (int L, int W, int x1, int x2, int y1, int y2)
{
  CutPoint c = { x1, x2, y1, y2, 0 };
  cutPoints[indexX[L]][indexY[W]] = c;
}

/******************************************************************
 ******************************************************************/

/**
 * Auxiliar function.
 *
 * Parameters:
 *
 * L - Length of the pallet.
 * W - Width of the pallet.
 * l - Length of the boxes.
 * w - Width of the boxes.
 * n - Maximum search depth.
 *
 * numBlocks - Number of blocks in the current division of the pallet.
 *
 * L_, W_ - Lenghts and widths for each partition of the pallet.
 *
 * z_lb   - Current lower bound for (L,W).
 * z_ub   - Upper bound for (L,W).
 *
 * x1, x2, y1, y2 - Points that determine the division of the pallet.
 *
 * Return:
 * - the number of (l,w)-boxes packed into (L,W) pallet, using the
 *   division determined by (x1,x2,y1,y2).
 */
int
solve (int L, int W, int l, int w, int n, int numBlocks, int *L_, int *W_,
       int *z_lb, int z_ub, int x1, int x2, int y1, int y2)
{

  /* z[1..5] stores the amount of boxes packed into partitions 1 to 5. */
  int z[6];

  /* Sum of the lower and upper bounds for the number of boxes that
   * can be packed into partitions 1 to 5. */
  int S_lb, S_ub;

  /* Indices of the subproblems in the indexing matrices. */
  int iX[6], iY[6];

  /* Lower and upper bounds for each partition. */
  int zi_ub[6], zi_lb[6];

  int i;

  /* Normalize each rectangle produced. */
  for (i = 1; i <= numBlocks; i++)
    {

      /* Normalize the size of the rectangle (Li,Wi). */
      L_[i] = normalize[L_[i]];
      W_[i] = normalize[W_[i]];

      /* We assume that Li >= Wi. */
      if (L_[i] < W_[i])
        {
          std::swap (L_[i], W_[i]);
        }

      /* Get the indices of each subproblem in the indexing matrices. */
      iX[i] = indexX[L_[i]];
      iY[i] = indexY[W_[i]];
    }

  /* If maximum level of the recursion was not reached. */
  if (n < N)
    {

      /* Store the sum of best packing estimations in the 5 partitions
       * until this moment. Initially, it receives the sum of the lower
       * bounds. */
      S_lb = 0;

      /* Sum of the upper bounds in the 5 partitions. */
      S_ub = 0;

      /* Compute the lower (zi_lb) and upper (zi_ub) bounds of each
       * partition. */
      for (i = 1; i <= numBlocks; i++)
        {
          /* Lower bound of (Li, Wi). */
          zi_lb[i] = lowerBound[iX[i]][iY[i]];
          S_lb += zi_lb[i];
          /* Upper bound of (Li, Wi). */
          zi_ub[i] = localUpperBound (iX[i], iY[i]);
          S_ub += zi_ub[i];
        }

      if (*z_lb < S_ub)
        {
          /* The current lower bound is less than the sum of the partitions
           * upper bounds. Then, there is a possibility of this division
           * improve the solution. */

          for (i = 1; i <= numBlocks; i++)
            {

#ifdef N_INFINITY
              if (solutionDepth[iX[i]][iY[i]] >= 0)
                {
                  z[i] = BD (L_[i], W_[i], l, w, n + 1);
                  lowerBound[iX[i]][iY[i]] = z[i];
                  solutionDepth[iX[i]][iY[i]] = -1;
                }
              else
                {
                  /* The problem was already solved. */
                  z[i] = lowerBound[iX[i]][iY[i]];
                }
#else
              if (solutionDepth[iX[i]][iY[i]] > n)
                {
                  /* Solve for the first time or give another chance for
                   * this problem. */
                  z[i] = BD (L_[i], W_[i], l, w, n + 1);
                  lowerBound[iX[i]][iY[i]] = z[i];

                  if (reachedLimit[iX[i]][iY[i]] == 0)
                    {
                      solutionDepth[iX[i]][iY[i]] = -1;
                    }
                  else
                    {
                      solutionDepth[iX[i]][iY[i]] = n;
                    }
                }
              else
                {
                  /* This problem was already solved. It gets the
                   * solution obtained previously. */
                  z[i] = lowerBound[iX[i]][iY[i]];
                }
#endif

              if (reachedLimit[iX[i]][iY[i]] == 1)
                {
                  reachedLimit[indexX[L]][indexY[W]] = 1;
                }

              /* Update lower and upper bounds for this partitioning. */
              S_lb = S_lb - zi_lb[i] + z[i];
              S_ub = S_ub - zi_ub[i] + z[i];

              /* If z_lb >= S_ub, we have, at least, a solution as good as
               * the one that can be find with this partitioning. So this
               * partitioning is discarded. */
              if (*z_lb >= S_ub)
                {
                  break;
                }

              /* If the sum of packings in the current partitions is better
               * than the previous, update the lower bound. */
              else if (S_lb > *z_lb)
                {
                  *z_lb = S_lb;
                  storeCutPoint (L, W, x1, x2, y1, y2);
                  if (*z_lb == z_ub)
                    {
                      /* An optimal solution was found. */
                      solutionDepth[indexX[L]][indexY[W]] = -1;
                      reachedLimit[indexX[L]][indexY[W]] = 0;
                      return 1;
                    }
                }
            }
        }
    } /* if n < N */

  /* The maximum depth of recursion was reached. Then, each partition
   * will not be solved recursivelly. Each one receives the best
   * packing obtained until this moment. */
  else
    {

      reachedLimit[indexX[L]][indexY[W]] = 1;
      S_lb = 0;

      /* Compute the lower bound of each partition and the sum is
       * stored in S_lb. */
      for (i = 1; i <= numBlocks; i++)
        {
          S_lb += lowerBound[iX[i]][iY[i]];
        }

      /* If the sum of the homogeneous packing in all current
       * partitions is better than the previous estimation for (L,W),
       * update the lower bound. */
      if (S_lb > *z_lb)
        {
          *z_lb = S_lb;
          storeCutPoint (L, W, x1, x2, y1, y2);
          if (*z_lb == z_ub)
            {
              /* An optimal solution was found. */
              solutionDepth[indexX[L]][indexY[W]] = -1;
              reachedLimit[indexX[L]][indexY[W]] = 0;
              return 1;
            }
        }
    }
  return 0;
}

/******************************************************************
 ******************************************************************/

/**
 * Guillotine and first order non-guillotine cuts recursive procedure.
 *
 * Parameters:
 * L - Length of the pallet.
 * W - Width of the pallet.
 * l - Length of the boxes.
 * w - Width of the boxes.
 * n - Maximum search depth.
 *
 * Return:
 * - the number of (l,w)-boxes packed into (L,W) pallet.
 */
int
BD (int L, int W, int l, int w, int n)
{

  /* Lower and upper bounds for the number of (l,w)-boxes that can be
   * packed into (L,W). */
  int z_lb, z_ub;

  /* We assume that L >= W. */
  if (W > L)
    {
      std::swap (L, W);
    }

  z_lb = lowerBound[indexX[L]][indexY[W]];
  z_ub = localUpperBound (indexX[L], indexY[W]);

  if (z_lb == 0 || z_lb == z_ub)
    {
      /* An optimal solution was found: no box fits into the pallet or
       * lower and upper bounds are equal. */
      solutionDepth[indexX[L]][indexY[W]] = -1;
      reachedLimit[indexX[L]][indexY[W]] = 0;
      return z_lb;
    }
  else
    {
      /* Points that determine the pallet division. */
      int x1, x2, y1, y2;

      /* Indices of x1, x2, y1 and y2 in the raster points arrays. */
      int index_x1, index_x2, index_y1, index_y2;

      /* Size of the generated partitions:
       * (L_[i], W_[i]) is the size of the partition i, for i = 1, ..., 5. */
      int L_[6], W_[6];

      /* Raster points sets X' and Y' for this problem. */
      Set rasterX, rasterY;

      int solved;

      /* Construct the raster points sets. */
      constructRasterPoints (L, W, &rasterX, &rasterY, normalSetX);

      reachedLimit[indexX[L]][indexY[W]] = 0;

      /*
       * Loop to generate the cut points (x1, x2, y1 and y2) considering
       * the following symmetries.
       *
       * - First order non-guillotine cuts:
       *   0 < x1 < x2 < L and 0 < y1 < y2 < W
       *   (x1 + x2 < L) or (x1 + x2 = L and y1 + y2 <= W)
       *
       * - Vertical guillotine cuts:
       *   0 < x1 = x2 <= L/2 and y1 = y2 = 0
       *
       * - Horizontal guillotine cuts:
       *   0 < y1 = y2 <= W/2 and x1 = x2 = 0
       *
       * Observation: In the loop of non-guillotine cuts it can appear less
       * than five partitions in the case of the normalization of a side
       * of a partition be zero.
       */

      /*#################################*
       * First order non-guillotine cuts *
       *#################################*/

      /*
           L_1     L_2
          ----------------
         |     |    2     |W_2
      W_1|  1  |          |
         |     |----------|
         |     | 3 |      |
         |---------|      |
         |         |  5   |W_5
      W_4|    4    |      |
         |         |      |
          ----------------
             L_4     L_5
      */

      for (index_x1 = 1;
           index_x1 < rasterX.size && rasterX.points[index_x1] <= L / 2;
           index_x1++)
        {

          x1 = rasterX.points[index_x1];

          for (index_x2 = index_x1 + 1;
               index_x2 < rasterX.size && rasterX.points[index_x2] + x1 <= L;
               index_x2++)
            {

              x2 = rasterX.points[index_x2];

              for (index_y1 = 1;
                   index_y1 < rasterY.size && rasterY.points[index_y1] < W;
                   index_y1++)
                {

                  y1 = rasterY.points[index_y1];

                  for (index_y2 = index_y1 + 1;
                       index_y2 < rasterY.size && rasterY.points[index_y2] < W;
                       index_y2++)
                    {

                      y2 = rasterY.points[index_y2];

                      /* Symmetry. When x1 + x2 = L, we can restrict y1 and y2
                       * to y1 + y2 <= W. */
                      if (x1 + x2 == L && y1 + y2 > W)
                        break;

                      /* The five partitions. */
                      L_[1] = x1;
                      W_[1] = W - y1;

                      L_[2] = L - x1;
                      W_[2] = W - y2;

                      L_[3] = x2 - x1;
                      W_[3] = y2 - y1;

                      L_[4] = x2;
                      W_[4] = y1;

                      L_[5] = L - x2;
                      W_[5] = y2;

                      solved = solve (L, W, l, w, n, 5, L_, W_, &z_lb, z_ub,
                                      x1, x2, y1, y2);

                      if (solved)
                        {
                          /* This problem was solved with optimality guarantee.
                           */
                          free (rasterX.points);
                          free (rasterY.points);
                          return z_lb;
                        }
                    } /* for y2 */
                }     /* for y1 */
            }         /* for x2 */
        }             /* for x1 */

      /*###########################*
       * Vertical guillotine cuts. *
       *###########################*/

      /*
         ----------------
        |     |          |
        |     |          |
        |     |          |
        |  1  |    2     |
        |     |          |
        |     |          |
        |     |          |
        |     |          |
         ----------------
      */

      for (index_x1 = 1;
           index_x1 < rasterX.size && rasterX.points[index_x1] <= L / 2;
           index_x1++)
        {

          x1 = x2 = rasterX.points[index_x1];
          y1 = y2 = 0;

          /* Partitions 1 and 2 generated by the vertical cut. */
          L_[1] = x1;
          W_[1] = W - y1;

          L_[2] = L - x1;
          W_[2] = W - y2;

          solved
              = solve (L, W, l, w, n, 2, L_, W_, &z_lb, z_ub, x1, x2, y1, y2);
          if (solved)
            {
              /* This problem was solved with optimality guarantee. */
              free (rasterX.points);
              free (rasterY.points);
              return z_lb;
            }
        }

      /*#############################*
       * Horizontal guillotine cuts. *
       *#############################*/

      /*
         ----------------
        |       2        |
        |                |
        |----------------|
        |                |
        |                |
        |       5        |
        |                |
        |                |
         ----------------
      */

      for (index_y1 = 1;
           index_y1 < rasterY.size && rasterY.points[index_y1] <= W / 2;
           index_y1++)
        {

          y1 = y2 = rasterY.points[index_y1];
          x1 = x2 = 0;

          /* Partitions 2 and 5 generated by the horizontal cut. */
          L_[1] = L - x1;
          W_[1] = W - y2;

          L_[2] = L - x2;
          W_[2] = y2;

          solved
              = solve (L, W, l, w, n, 2, L_, W_, &z_lb, z_ub, x1, x2, y1, y2);
          if (solved)
            {
              /* This problem was solved with optimality guarantee. */
              free (rasterX.points);
              free (rasterY.points);
              return z_lb;
            }
        }

      free (rasterX.points);
      free (rasterY.points);
      return z_lb;
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Compute the Barnes's upper bound [3].
 *
 * [3] F. W. Barnes. Packing the maximum number of m x n tiles in a
 *     large p x q rectangle. Discrete Mathematics, volume 26,
 *     pages 93-100, 1979.
 *
 * Parameters:
 * (L, W) Pallet dimensions.
 * (l, w) Dimensions of the boxes to be packed.
 *
 * Return:
 * - the computed Barnes's bound.
 */
int
barnesBound (int L, int W, int l, int w)
{
  int r, s, D;
  int minWaste = (L * W) % (l * w);

  /* (l,1)-boxes packing. */
  r = L % l;
  s = W % l;
  int A = std::min (r * s, (l - r) * (l - s));

  /* (1,w)-boxes packing. */
  r = L % w;
  s = W % w;
  int B = std::min (r * s, (w - r) * (w - s));

  /* Best unitary tile packing. */
  int maxAB = std::max (A, B);

  if (minWaste >= maxAB % (l * w))
    {
      /* Wasted area. */
      D = (maxAB / (l * w)) * (l * w) + minWaste;
    }
  else
    {
      /* Wasted area. */
      D = (maxAB / (l * w) + 1) * (l * w) + minWaste;
    }

  return (L * W - D) / (l * w);
}

/******************************************************************
 ******************************************************************/

void
initialize (int L, int W, int l, int w)
{

  /* Normalization of L and W. */
  int L_n, W_n;

  int i, j;

  Set rasterX, rasterY;

  /* Construct the conic combination set of l and w. */
  normalSetX = newSet (L + 2);
  constructConicCombinations (L, l, w, &normalSetX);

  /* Compute the values of L* and W*.
   * normalize[i] = max {x in X | x <= i} */
  normalize = new int[L + 1];
  i = 0;
  for (j = 0; j <= L; j++)
    {
      for (; i < normalSetX.size && normalSetX.points[i] <= j; i++)
        ;
      normalize[j] = normalSetX.points[i - 1];
    }

  /* Normalize (L, W). */
  L_n = normalize[L];
  W_n = normalize[W];

  constructRasterPoints (L, W, &rasterX, &rasterY, normalSetX);

  free (normalSetX.points);
  normalSetX = newSet (L_n + 2);
  int k = 0;
  i = 0;
  j = 0;
  while (i < rasterX.size && rasterX.points[i] <= L_n && j < rasterY.size
         && rasterY.points[j] <= W_n)
    {

      if (rasterX.points[i] == rasterY.points[j])
        {
          normalSetX.points[k++] = rasterX.points[i++];
          normalSetX.size++;
          j++;
        }
      else if (rasterX.points[i] < rasterY.points[j])
        {
          normalSetX.points[k++] = rasterX.points[i++];
          normalSetX.size++;
        }
      else
        {
          normalSetX.points[k++] = rasterY.points[j++];
          normalSetX.size++;
        }
    }
  while (i < rasterX.size && rasterX.points[i] <= L_n)
    {
      if (rasterX.points[i] > normalSetX.points[k - 1])
        {
          normalSetX.points[k++] = rasterX.points[i];
          normalSetX.size++;
        }
      i++;
    }
  if (k > 0 && normalSetX.points[k - 1] < L_n)
    {
      normalSetX.points[k++] = L_n;
      normalSetX.size++;
    }
  normalSetX.points[k] = L_n + 1;
  normalSetX.size++;

  free (rasterX.points);
  free (rasterY.points);

  /* Construct the array of indices. */
  indexX = new int[L_n + 2];
  indexY = new int[W_n + 2];

  for (i = 0; i < normalSetX.size; i++)
    {
      indexX[normalSetX.points[i]] = i;
    }
  int ySize = 0;
  for (i = 0; i < normalSetX.size; i++)
    {
      if (normalSetX.points[i] > W_n)
        {
          break;
        }
      ySize++;
      indexY[normalSetX.points[i]] = i;
    }

  solutionDepth = new int *[normalSetX.size];
  upperBound = new int *[normalSetX.size];
  lowerBound = new int *[normalSetX.size];
  reachedLimit = new int *[normalSetX.size];
  cutPoints = new CutPoint *[normalSetX.size];

  for (i = 0; i < normalSetX.size; i++)
    {
      solutionDepth[i] = new int[ySize];
      upperBound[i] = new int[ySize];
      lowerBound[i] = new int[ySize];
      cutPoints[i] = new CutPoint[ySize];
      reachedLimit[i] = new int[ySize];
    }

  for (i = 0; i < normalSetX.size; i++)
    {

      int x = normalSetX.points[i];

      for (j = 0; j < ySize; j++)
        {

          int y = normalSetX.points[j];

          solutionDepth[i][j] = N;
          reachedLimit[i][j] = 1;
          upperBound[i][j] = barnesBound (x, y, l, w);
          lowerBound[i][j] = lowerBound (x, y, l, w);
          cutPoints[i][j].homogeneous = 1;
        }
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Parameters:
 * L     - Length of the pallet.
 * W     - Width of the pallet.
 * l     - Length of the boxes.
 * w     - Width of the boxes.
 * N_max - Maximum depth.
 */
int
solve_BD (int L, int W, int l, int w, int N_max)
{
  int L_n, W_n;

  N = N_max;
  if (N <= 0)
    {
      N = INFINITY_;
    }

  /* We assume that L >= W. */
  if (W > L)
    {
      std::swap (L, W);
    }

  initialize (L, W, l, w);

  /* Normalize (L, W). */
  L_n = normalize[L];
  W_n = normalize[W];

  int solution = BD (L_n, W_n, l, w, N);

  /* remove this stupid check */
  // if (solution != upperBound[indexX[L_n]][indexY[W_n]] && N != 1)
  //   {
  //     /* If the solution found is not optimal (or, at least, if it is not
  //      * possible to prove its optimality), solve from the first level. */
  //     printf ("solving depth 1: %d \n", solution);
  //     solution = BD (L_n, W_n, l, w, 1);
  //   }

  lowerBound[indexX[L_n]][indexY[W_n]] = solution;

  for (int i = 0; i < normalSetX.size; i++)
    {
      free (reachedLimit[i]);
      free (solutionDepth[i]);
    }

  delete[] reachedLimit;
  delete[] solutionDepth;

  return solution;
}
