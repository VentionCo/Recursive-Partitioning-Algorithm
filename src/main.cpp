/* The Recursive Partitioning Algorithm is designed to solve the
 * manufacturer's pallet loading problem. This problem consists in
 * arranging, without overlapping, identical rectangular pieces with
 * length l and width w in a rectangular pallet with length L and
 * width W. The pieces must be placed orthogonally and only 90-degree
 * rotations are allowed. The objective is to find a pattern with the
 * maximum number of pieces packed.
 *
 * Copyright (C) 2007-2016 Rafael Durbano Lobato
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact information:
 *
 * Rafael Durbano Lobato <lobato@ime.usp.br>
 * http://www.ime.usp.br/~lobato/
 */

/* Implementation of the "Recursive Partitioning Approach" described
 * in [1] for solving the Manufacturer's Pallet Loading Problem.
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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>

#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/times.h>

#include <iostream>

#include "bd.h"
#include "draw.h"
#include "graphics.h"
#include "sets.h"
#include "util.h"

#include <emscripten.h>

// If this is an Emscripten (WebAssembly) build then...
#ifdef __EMSCRIPTEN__
  #include <emscripten.h>
#endif

/******************************************************************
 ******************************************************************/

int solve (int L, int *q);

/* Lower and upper bounds of each rectangular subproblem. */
int **lowerBound, **upperBound;

/* Arrays of indices for indexing the matrices that store information
 * about rectangular subproblems (L,W), where (L,W) belongs to X' x Y'
 * and X' and Y' are the raster points sets associated to (L,l,w) and
 * (W,l,w), respectively. */
int *indexX, *indexY;

/* Set of integer conic combinations of l and w:
 * X = {x | x = rl + sw, with r,w in Z and r,w >= 0} */
Set normalSetX;

/* Array that stores the normalized values of each integer between 0 and
 * L (dimension of the problem):
 * normalize[x] = max {r in X' | r <= x} */
int *normalize;

/* Store the solutions of the subproblems. */
std::map<int, int> *solutionMap;
int *solution;

/* Store the division points in the rectangular and in the L-shaped
 * pieces associated to the solutions found. */
std::map<int, int> *divisionPointMap;
int *divisionPoint;

/* Dimensions of the boxes to be packed. */
int l, w;

/* Type of the structure used to store the solutions. */
int memory_type;

/* Store the points that determine the divisions of the rectangles. */
CutPoint **cutPoints;

int *indexRasterX, *indexRasterY;
int numRasterX, numRasterY;

/******************************************************************
 ******************************************************************/

inline int
roundToNearest (double a)
{
  return (int)floor (a + 0.5);
}

/******************************************************************
 ******************************************************************/

/**
 * Calculate the upper bound of a given rectangle R(x,y) (degenerated L).
 *
 * Parameters:
 * x - Length of the rectangle.
 *
 * y - Width of the rectangle.
 *
 * Return:
 * The computed upper bound.
 */
inline int
R_UpperBound (int x, int y)
{
  /* A(R) / lw */
  x = normalize[x];
  y = normalize[y];
  return upperBound[indexX[x]][indexY[y]];
}

/******************************************************************
 ******************************************************************/

/**
 * Calculate the upper bound of a given L.
 *
 * Parameters:
 * q - The L-piece.
 *
 * Return:
 * The computed upper bound for this L-piece.
 */
inline int
L_UpperBound (int *q)
{
  /* Area(L) / lw */
  return (q[0] * q[1] - (q[0] - q[2]) * (q[1] - q[3])) / (l * w);
}

/******************************************************************
 ******************************************************************/

/**
 * Calculate the lower bound of a given rectangle R(x,y) (degenerated L).
 *
 * Parameters:
 * x - Length of the rectangle.
 *
 * y - Width of the rectangle.
 *
 * Return:
 * The computed lower bound.
 */
inline int
R_LowerBound (int x, int y)
{
  x = normalize[x];
  y = normalize[y];
  return lowerBound[indexX[x]][indexY[y]];
}

/******************************************************************
 ******************************************************************/

/**
 * Calculate the lower bound of a given L. It divides the L in two
 * rectangles and calculates their lower bounds to compose the lower
 * bound of the L-piece.
 *
 * +-----+              +-----+              +-----+
 * |     |              |     |              |     |
 * |     |              |     |              |     |
 * |     +----+   -->   +-----+----+    or   |     +----+
 * |          |         |          |         |     |    |
 * |          |         |          |         |     |    |
 * +----------+         +----------+         +-----+----+
 *                           (a)                  (b)
 *
 * Parameters:
 * q - The L-piece.
 *
 * Return:
 * The computed lower bound.
 */
inline int
L_LowerBound (int *q, bool *horizontalCut)
{

  int a = lowerBound[indexX[normalize[q[2]]]][indexY[normalize[q[1] - q[3]]]]
          + lowerBound[indexX[normalize[q[0]]]][indexY[normalize[q[3]]]];

  int b
      = lowerBound[indexX[normalize[q[2]]]][indexY[normalize[q[1]]]]
        + lowerBound[indexX[normalize[q[0] - q[2]]]][indexY[normalize[q[3]]]];

  if (a > b)
    {
      *horizontalCut = true;
      return a;
    }
  else
    {
      *horizontalCut = false;
      return b;
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Divide an L-piece in two new L-pieces, according to the specified
 * subdivision, and normalize the two ones.
 *
 * Parameters:
 * i  - Point that determines the division int the L-piece.
 *
 * q  - The L-piece to be divided.
 *
 * q1 - It will store a new L-piece.
 *
 * q2 - It will store the other new L-piece.
 *
 * standardPosition - Pointer to the function that will divide the L-piece.
 */
void
divide (int *i, int *q, int *q1, int *q2,
        void (*standardPosition) (int *, int *, int *, int *))
{

  /* Divide the L-piece in two new ones. */
  (*standardPosition) (i, q, q1, q2);

  /* Normalize the new L-pieces. */
  normalizePiece (q1);
  normalizePiece (q2);
}

/******************************************************************
 ******************************************************************/

/**
 * Return the solution of the L-piece related to the index L.
 *
 * Parameters:
 * L   - Index of the L-piece.
 *
 * key - Key for this L-piece.
 *
 * Return:
 * The current solution of the specified L-piece.
 */
inline int
getSolution (int L, int key)
{
  if (memory_type == MEM_TYPE_4)
    {
      return solution[L] & nRet;
    }
  else
    {
      return solutionMap[L][key] & nRet;
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Return the solution of the L-piece related to the index L.
 *
 * Parameters:
 * L - Index of the L-piece.
 *
 * q - The L-piece.
 *
 * Returns:
 * The current solution of the specified L-piece.
 */
inline int
getSolution (int L, int *q)
{
  if (memory_type == MEM_TYPE_4)
    {
      return solution[L] & nRet;
    }
  else
    {
      int key = getKey (q[0], q[1], q[2], q[3], memory_type);
      return solutionMap[L][key] & nRet;
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Return the solution of the L-piece related to the index L.
 *
 * Parameters:
 * L   - Index of the L-piece.
 *
 * q   - The L-piece.
 *
 * key - Key for this L-piece.
 *
 * Return:
 * The current solution of the specified L-piece.
 */
inline int
getSolution (int L, int *q, int *key)
{
  if (memory_type == MEM_TYPE_4)
    {
      return solution[L];
    }
  else
    {
      *key = getKey (q[0], q[1], q[2], q[3], memory_type);
      return solutionMap[L][*key];
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Store the solution of an L-piece.
 *
 * Parameters:
 * L         - Index of the L-piece.
 *
 * key       - Key for this L-piece.
 *
 * LSolution - Solution to be stored.
 *
 */
inline void
storeSolution (int L, int key, int LSolution)
{
  if (memory_type == MEM_TYPE_4)
    {
      solution[L] = LSolution;
    }
  else
    {
      solutionMap[L][key] = LSolution;
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Store the point where the division in the L-piece was made.
 *
 * Parameters:
 * L     - Index of the L-piece.
 *
 * key   - Key for this L-piece.
 *
 * point - Representation of the point where the division was made.
 *
 */
inline void
storeDivisionPoint (int L, int key, int point)
{
  if (memory_type == MEM_TYPE_4)
    {
      divisionPoint[L] = point;
    }
  else
    {
      divisionPointMap[L][key] = point;
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Verify whether the solution is optimal.
 *
 * Parameters:
 * L          - Index of the L-piece.
 *
 * key        - Key for this L-piece.
 *
 * upperBound - Upper bound for the L-piece.
 *
 * Return:
 * Return whether the current solution for the L-piece is optimal.
 */
inline bool
optimal (int L, int key, int upperBound)
{
  int Lsolution = getSolution (L, key);
  if ((Lsolution & nRet) == upperBound)
    {
      return true;
    }
  return false;
}

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-piece in every possible way, according to the specified
 * subdivision B.
 *
 * Parameters:
 * L                - Index of the L-piece.
 *
 * q                - The L-piece.
 *
 * constraints      - Constraints that determine the interval of x' and y'.
 *
 * B                - The subdivision.
 *
 * standardPosition - Pointer to the function that divides the L-piece
 *                    according to the subdivision B.
 *
 * X                - Set of raster points.
 *
 * startX           - Index to start the divisions on the set X.
 *
 * Y                - Set of raster points.
 *
 * startY           - Index to start the divisions on the set Y.
 */
int
divideL (int L, int *q, int *constraints, int B,
         void (*standardPosition) (int *, int *, int *, int *), Set X,
         int startX, Set Y, int startY)
{

  /* i_k[0] <- x'
   * i_k[1] <- y'
   */
  int i_k[2];
  int i_x, i_y;
  int q1[4], q2[4];

  int key = 0;
  int LSolution = getSolution (L, q, &key);
  int upperBound = L_UpperBound (q);

  for (i_x = startX; i_x < X.size; i_x++)
    {

      i_k[0] = X.points[i_x];
      if (i_k[0] > constraints[1])
        {
          break;
        }

      for (i_y = startY; i_y < Y.size; i_y++)
        {

          i_k[1] = Y.points[i_y];
          if (i_k[1] > constraints[3])
            {
              break;
            }

          divide (i_k, q, q1, q2, (*standardPosition));
          if (q1[0] < 0 || q2[0] < 0)
            {
              continue;
            }

          if (L_UpperBound (q1) + L_UpperBound (q2) > (LSolution & nRet))
            {
              /* It is possible that this division gets a better solution. */
              int L1 = LIndex (q1[0], q1[1], q1[2], q1[3], memory_type);
              int L2 = LIndex (q2[0], q2[1], q2[2], q2[3], memory_type);
              int L1Solution = solve (L1, q1);
              int L2Solution = solve (L2, q2);

              if ((L1Solution & nRet) + (L2Solution & nRet)
                  > (LSolution & nRet))
                {
                  /* A better solution was found. */
                  LSolution = ((L1Solution & nRet) + (L2Solution & nRet))
                              | (B << descSol);
                  storeSolution (L, key, LSolution);
                  storeDivisionPoint (L, key,
                                      i_k[0] | (i_k[1] << descPtoDiv2));
                  if ((LSolution & nRet) == upperBound)
                    {
                      return LSolution;
                    }
                }
            }
        }
    }
  return LSolution;
}

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-piece in every possible way, according to the B6
 * subdivision.
 *
 * +-------------+--------+
 * |             |        |
 * |   (x',y')   |   L2   |
 * |      o------o        |
 * |      |  (x'',y')     |
 * |  L1  |               |
 * |      |               |
 * +------+---------------+
 *
 * Parameters:
 * L - Index of the L-piece.
 *
 * q - The L-piece.
 *
 * X - Set of raster points.
 *
 * Y - Set of raster points.
 */
int
divideB6 (int L, int *q, Set X, Set Y)
{

  /* i_k[0] <- x'
   * i_k[1] <- y'
   * i_k[2] <- x''
   */
  int i_k[3];
  int q1[4], q2[4];
  int key = 0;
  int LSolution = getSolution (L, q, &key);
  int upperBound = R_UpperBound (q[0], q[1]);

  int i = 0;
  for (i_k[0] = X.points[i]; i < X.size; i++)
    {

      i_k[0] = X.points[i];

      int j = i;
      for (i_k[2] = X.points[j]; j < X.size; j++)
        {

          i_k[2] = X.points[j];
          if (i_k[0] == 0 && i_k[2] == 0)
            {
              continue;
            }

          int k = 0;
          for (i_k[1] = Y.points[k]; k < Y.size; k++)
            {

              i_k[1] = Y.points[k];
              divide (i_k, q, q1, q2, standardPositionB6);
              if (q1[0] < 0 || q2[0] < 0)
                {
                  continue;
                }

              if (L_UpperBound (q1) + L_UpperBound (q2) > (LSolution & nRet))
                {
                  /* It is possible that this division gets a better solution.
                   */
                  int L1 = LIndex (q1[0], q1[1], q1[2], q1[3], memory_type);
                  int L2 = LIndex (q2[0], q2[1], q2[2], q2[3], memory_type);
                  int L1Solution = solve (L1, q1);
                  int L2Solution = solve (L2, q2);

                  if ((L1Solution & nRet) + (L2Solution & nRet)
                      > (LSolution & nRet))
                    {
                      /* A better solution was found. */
                      LSolution = ((L1Solution & nRet) + (L2Solution & nRet))
                                  | (B6 << descSol);
                      storeSolution (L, key, LSolution);
                      storeDivisionPoint (L, key,
                                          i_k[0] | (i_k[1] << descPtoDiv2)
                                              | (i_k[2] << descPtoDiv3));

                      if ((LSolution & nRet) == upperBound)
                        {
                          return LSolution;
                        }
                    }
                }
            }
        }
    }
  return LSolution;
}

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-piece in every possible way, according to the B7
 * subdivision.
 *
 * +-------------+
 * |             |
 * |   (x',y'')  |
 * |      o------+
 * |      |      |
 * |  L1  |  L2  |
 * |      |      |
 * +------o      |
 * |   (x',y')   |
 * |             |
 * |             |
 * +-------------+
 *
 * Parameters:
 * L - Index of the L-piece.
 *
 * q - The L-piece.
 *
 * X - Set of raster points.
 *
 * Y - Set of raster points.
 */
int
divideB7 (int L, int *q, Set X, Set Y)
{

  /* i_k[0] <- x'
   * i_k[1] <- y'
   * i_k[2] <- y''
   */
  int i_k[3];
  int q1[4], q2[4];
  int key = 0;
  int LSolution = getSolution (L, q, &key);
  int upperBound = R_UpperBound (q[0], q[1]);

  int j = 0;
  for (i_k[1] = Y.points[j]; j < Y.size; j++)
    {

      i_k[1] = Y.points[j];

      int k = j;
      for (i_k[2] = Y.points[k]; k < Y.size; k++)
        {

          i_k[2] = Y.points[k];
          if (i_k[1] == 0 && i_k[2] == 0)
            {
              continue;
            }

          int i = 0;
          for (i_k[0] = X.points[i]; i < X.size; i++)
            {

              i_k[0] = X.points[i];
              divide (i_k, q, q1, q2, standardPositionB7);
              if (q1[0] < 0 || q2[0] < 0)
                {
                  continue;
                }

              if (L_UpperBound (q1) + L_UpperBound (q2) > (LSolution & nRet))
                {
                  /* It is possible that this division gets a better solution.
                   */
                  int L1 = LIndex (q1[0], q1[1], q1[2], q1[3], memory_type);
                  int L2 = LIndex (q2[0], q2[1], q2[2], q2[3], memory_type);
                  int L1Solution = solve (L1, q1);
                  int L2Solution = solve (L2, q2);

                  if (((L1Solution & nRet) + (L2Solution & nRet))
                      > (LSolution & nRet))
                    {
                      /* A better solution was found. */
                      LSolution = ((L1Solution & nRet) + (L2Solution & nRet))
                                  | (B7 << descSol);
                      storeSolution (L, key, LSolution);
                      storeDivisionPoint (L, key,
                                          i_k[0] | (i_k[1] << descPtoDiv2)
                                              | (i_k[2] << descPtoDiv3));

                      if ((LSolution & nRet) == upperBound)
                        {
                          return LSolution;
                        }
                    }
                }
            }
        }
    }
  return LSolution;
}

/******************************************************************
 ******************************************************************/

/**
 * Solve the problem of packing rectangular (l,w)-boxes into the
 * specified L-shaped piece.
 *
 * Parameters:
 * L - Index of the L-piece.
 *
 * q - The L-piece. q = {X, Y, x, y}.
 */
int
solve (int L, int *q)
{

  int key = 0;
  if (memory_type == MEM_TYPE_4)
    {
      if (solution[L] != -1)
        {
          /* This problem has already been solved. */
          return solution[L];
        }
    }
  else
    {
      key = getKey (q[0], q[1], q[2], q[3], memory_type);
      if (solutionMap[L].count (key) > 0)
        {
          /* This problem has already been solved. */
          return solutionMap[L][key];
        }
    }

  if (q[0] != q[2])
    {
      bool horizontalCut;
      int lowerBound = L_LowerBound (q, &horizontalCut);
      int upperBound = L_UpperBound (q);
      int LSolution = lowerBound | (B1 << descSol);

      if (horizontalCut)
        storeDivisionPoint (L, key, 0 | (q[3] << descPtoDiv2));
      else
        storeDivisionPoint (L, key, q[2] | (0 << descPtoDiv2));

      storeSolution (L, key, LSolution);

      /* Try to solve this problem with homogeneous packing (or other
       * better solution already computed). */
      if ((LSolution & nRet) != upperBound)
        {
          /* It was not possible to solve this problem with homogeneous
           * packing. */
          int constraints[4];
          int startX = 0;
          int startY = 0;

          /* Construct the raster points sets X and Y. */
          Set X, Y;
          constructRasterPoints (q[0], q[1], &X, &Y, normalSetX);
          for (startX = 0; X.points[startX] < q[2]; startX++)
            ;
          for (startY = 0; Y.points[startY] < q[3]; startY++)
            ;

          /***********************************
           * 0 <= x' <= x  and  0 <= y' <= y *
           ***********************************/
          constraints[0] = 0;
          constraints[1] = q[2];
          constraints[2] = 0;
          constraints[3] = q[3];

          /* B1 subdivision.
           *
           * +------------+
           * |            |
           * |            |(x,y)
           * |      +-----o-----+
           * |  L1  |           |
           * |      |     L2    |
           * +------o           |
           * |   (x',y')        |
           * |                  |
           * +------------------+
           */
          LSolution = divideL (L, q, constraints, B1, standardPositionB1, X, 0,
                               Y, 0);
          if ((LSolution & nRet) == upperBound)
            {
              free (X.points);
              free (Y.points);
              return LSolution;
            }

          /* B3 subdivision.
           *
           * +------+-----+
           * |      |     |
           * |      |     |(x,y)
           * |      | L2  o-----+
           * |      |           |
           * |  L1  |           |
           * |      o-----------+
           * |   (x',y')        |
           * |                  |
           * +------------------+
           */
          LSolution = divideL (L, q, constraints, B3, standardPositionB3, X, 0,
                               Y, 0);
          if ((LSolution & nRet) == upperBound)
            {
              free (X.points);
              free (Y.points);
              return LSolution;
            }

          /* B5 subdivision.
           *
           * +------------+
           * |            |
           * |     L1     |(x,y)
           * |            o-----+
           * |   (x',y')  |     |
           * |      o-----+     |
           * |      |           |
           * |      |     L2    |
           * |      |           |
           * +------+-----------+
           */
          LSolution = divideL (L, q, constraints, B5, standardPositionB5, X, 0,
                               Y, 0);
          if ((LSolution & nRet) == upperBound)
            {
              free (X.points);
              free (Y.points);
              return LSolution;
            }

          /***********************************
           * 0 <= x' <= x  and  y <= y' <= Y *
           ***********************************/
          constraints[0] = 0;
          constraints[1] = q[2];
          constraints[2] = q[3];
          constraints[3] = Y.points[Y.size - 1];

          /* B2 subdivision.
           *
           * +------------+
           * |            |
           * |   (x',y')  |
           * +------o     |
           * |      | L1  |
           * |      |     |(x,y)
           * |      +-----o-----+
           * |  L2              |
           * |                  |
           * +------------------+
           */
          LSolution = divideL (L, q, constraints, B2, standardPositionB2, X, 0,
                               Y, startY);
          if ((LSolution & nRet) == upperBound)
            {
              free (X.points);
              free (Y.points);
              return LSolution;
            }

          /* B8 subdivision.
           *
           * +------------+
           * |            |
           * |   (x',y')  |
           * |      o-----+
           * |      |     |
           * |  L1  |     |(x,y)
           * |      |     o-----+
           * |      |  L2       |
           * |      |           |
           * +------+-----------+
           */
          LSolution = divideL (L, q, constraints, B8, standardPositionB8, X, 0,
                               Y, startY);
          if ((LSolution & nRet) == upperBound)
            {
              free (X.points);
              free (Y.points);
              return LSolution;
            }

          /***********************************
           * x <= x' <= X  and  0 <= y' <= y *
           ***********************************/
          constraints[0] = q[2];
          constraints[1] = X.points[X.size - 1];
          constraints[2] = 0;
          constraints[3] = q[3];

          /* B4 subdivision.
           *
           * +------+
           * |      |
           * |      |(x,y)
           * |      o-----------+
           * |  L1  |           |
           * |      |  (x',y')  |
           * |      +-----o     |
           * |            | L2  |
           * |            |     |
           * +------------+-----+
           */
          LSolution = divideL (L, q, constraints, B4, standardPositionB4, X,
                               startX, Y, 0);
          if ((LSolution & nRet) == upperBound)
            {
              free (X.points);
              free (Y.points);
              return LSolution;
            }

          /* B9 subdivision.
           *
           * +---------+
           * |         |
           * |         |(x,y)
           * |   L1    o---+----+
           * |             |    |
           * |             |    |
           * +-------------o    |
           * |          (x',y') |
           * |     L2           |
           * |                  |
           * +------------------+
           */
          LSolution = divideL (L, q, constraints, B9, standardPositionB9, X,
                               startX, Y, 0);
          free (X.points);
          free (Y.points);
        }
      return LSolution;
    } /* if q[0] != q[2] */
  else
    {
      /* Degenerated L (a rectangle) */

      int LSolution = R_LowerBound (q[0], q[1]) | (HOMOGENEOUS << descSol);
      int upperBound = R_UpperBound (q[0], q[1]);
      storeSolution (L, key, LSolution);

      /* Verify whether it could not be solved with homogeneous packing. */
      if ((LSolution & nRet) != upperBound)
        {

          /* Construct the raster points sets X and Y. */
          Set X, Y;
          constructRasterPoints (q[0], q[1], &X, &Y, normalSetX);

          /* Try the subdivisions B6 and B7. */

          /* B6 subdivision.
           *
           * +-------------+--------+
           * |             |        |
           * |   (x',y')   |   L2   |
           * |      o------o        |
           * |      |  (x'',y')     |
           * |  L1  |               |
           * |      |               |
           * +------+---------------+
           */
          LSolution = divideB6 (L, q, X, Y);
          if ((LSolution & nRet) == upperBound)
            {
              free (X.points);
              free (Y.points);

              /* Update the lower bound for this rectangular piece. */
              lowerBound[indexX[q[0]]][indexY[q[1]]] = LSolution & nRet;

              return LSolution;
            }

          /* B7 subdivision.
           *
           * +-------------+
           * |             |
           * |   (x',y'')  |
           * |      o------+
           * |      |      |
           * |  L1  |  L2  |
           * |      |      |
           * +------o      |
           * |   (x',y')   |
           * |             |
           * |             |
           * +-------------+
           */
          LSolution = divideB7 (L, q, X, Y);

          free (X.points);
          free (Y.points);

          /* Update the lower bound for this rectangular piece. */
          lowerBound[indexX[q[0]]][indexY[q[1]]] = LSolution & nRet;
        }
      return LSolution;
    }
}

/******************************************************************
 ******************************************************************/

void
makeIndices (int L, int W)
{

  Set X, Y, raster;
  constructRasterPoints (L, W, &X, &Y, normalSetX);

  int j = 0;
  int k = 0;
  int i = 0;

  raster = newSet (L + 2);
  while (i < X.size && X.points[i] <= L && j < Y.size && Y.points[j] <= W)
    {

      if (X.points[i] == Y.points[j])
        {
          raster.points[k++] = X.points[i++];
          raster.size++;
          j++;
        }
      else if (X.points[i] < Y.points[j])
        {
          raster.points[k++] = X.points[i++];
          raster.size++;
        }
      else
        {
          raster.points[k++] = Y.points[j++];
          raster.size++;
        }
    }
  while (i < X.size && X.points[i] <= L)
    {
      if (X.points[i] > raster.points[k - 1])
        {
          raster.points[k++] = X.points[i];
          raster.size++;
        }
      i++;
    }
  if (k > 0 && raster.points[k - 1] < L)
    {
      raster.points[k++] = L;
      raster.size++;
    }
  raster.points[k] = L + 1;
  raster.size++;

  try
    {
      indexRasterX = new int[L + 2];
      indexRasterY = new int[W + 2];
    }
  catch (std::exception &e)
    {
      std::cout << "Error allocating memory." << std::endl;
      exit (0);
    }

  j = 0;
  numRasterX = 0;
  for (int i = 0; i <= L; i++)
    {

      if (raster.points[j] == i)
        {
          indexRasterX[i] = numRasterX++;
          j++;
        }
      else
        {
          indexRasterX[i] = indexRasterX[i - 1];
        }
    }
  indexRasterX[L + 1] = indexRasterX[L] + 1;

  j = 0;
  numRasterY = 0;
  for (int i = 0; i <= W; i++)
    {
      if (raster.points[j] == i)
        {
          indexRasterY[i] = numRasterY++;
          j++;
        }
      else
        {
          indexRasterY[i] = indexRasterY[i - 1];
        }
    }
  indexRasterY[W + 1] = indexRasterY[W] + 1;

  free (raster.points);
  free (X.points);
  free (Y.points);
}

/******************************************************************
 ******************************************************************/

void
freeMemory ()
{
  if (memory_type == MEM_TYPE_4)
    {
      delete[] solution;
      delete[] divisionPoint;
    }
  else
    {
      delete[] solutionMap;
      delete[] divisionPointMap;
    }
  delete[] indexRasterX;
  delete[] indexRasterY;
}

/******************************************************************
 ******************************************************************/

bool
tryAllocateMemory (int size)
{
  try
    {
      solutionMap = new std::map<int, int>[size];
    }
  catch (std::exception &e)
    {
      if (size == 0)
        {
          std::cout << "Error allocating memory." << std::endl;
          exit (0);
        }
      return false;
    }
  try
    {
      divisionPointMap = new std::map<int, int>[size];
    }
  catch (std::exception &e)
    {
      delete[] solutionMap;
      delete[] divisionPointMap;
      if (size == 0)
        {
          std::cout << "Error allocating memory." << std::endl;
          exit (0);
        }
      return false;
    }
  return true;
}

/******************************************************************
 ******************************************************************/

void
allocateMemory ()
{
  memory_type = MEM_TYPE_4;

  int nL = roundToNearest (
      (pow ((double)numRasterX, ceil ((double)memory_type / 2.0))
       * pow ((double)numRasterY, floor ((double)memory_type / 2.0))));

  memory_type--;

  if (nL >= 0)
    {
      try
        {
          solution = new int[nL];
          try
            {
              divisionPoint = new int[nL];
              for (int i = 0; i < nL; i++)
                solution[i] = -1;
            }
          catch (std::exception &e)
            {

              delete[] solution;
              do
                {
                  nL = roundToNearest (
                      (pow ((double)numRasterX,
                            ceil ((double)memory_type / 2.0))
                       * pow ((double)numRasterY,
                              floor ((double)memory_type / 2.0))));

                  memory_type--;

                  if (nL >= 0 && tryAllocateMemory (nL))
                    {
                      break;
                    }
                }
              while (memory_type >= 0);
            }
        }
      catch (std::exception &e)
        {
          do
            {
              nL = roundToNearest (
                  (pow ((double)numRasterX, ceil ((double)memory_type / 2.0))
                   * pow ((double)numRasterY,
                          floor ((double)memory_type / 2.0))));

              memory_type--;

              if (nL >= 0 && tryAllocateMemory (nL))
                {
                  break;
                }
            }
          while (memory_type >= 0);
        }
    }
  else
    {
      do
        {
          nL = roundToNearest (
              (pow ((double)numRasterX, ceil ((double)memory_type / 2.0))
               * pow ((double)numRasterY, floor ((double)memory_type / 2.0))));

          memory_type--;

          if (nL >= 0 && tryAllocateMemory (nL))
            {
              break;
            }
        }
      while (memory_type >= 0);
    }
  memory_type++;
}

/******************************************************************
 ******************************************************************/

#ifdef __cplusplus
extern "C" { // So that the C++ compiler does not rename our function names
#endif

#ifdef __EMSCRIPTEN__
  EMSCRIPTEN_KEEPALIVE
#endif
  const char* pack(int inL, int inW, int inl, int inw) {
    printf("Beginning pack sequence\n");
    int L, W;
    int q[4];
    int BD_solution;
    int L_n, W_n;
    bool swap = false;

    // int INDEX, L_solution;
    L = inL;
    W = inW;
    l = inl;
    w = inw;

    if (L <= 0 || W <= 0 || l <= 0 || w <= 0) {
      printf ("Invalid dimensions\n");
      return NULL;
    }

    if (L < W) {
      std::swap (L, W);
      swap = true;
    }

    memory_type = 5;

    /* Try to solve the problem with Algorithm 1. */
    BD_solution = solve_BD (L, W, l, w, 0);

    L_n = normalize[L];
    W_n = normalize[W];

    q[0] = q[2] = L_n;
    q[1] = q[3] = W_n;

    auto result = draw (L, W, 0, q, BD_solution, false, l, w, swap);

    for (int i = 0; i < normalSetX.size; i++) {
      free (lowerBound[i]);
      free (upperBound[i]);
      free (cutPoints[i]);
    }

    delete[] lowerBound;
    delete[] upperBound;
    delete[] cutPoints;
    delete[] indexX;
    delete[] indexY;
    delete[] normalize;
    delete[] normalSetX.points;

    return result.c_str();
  }

#ifdef __cplusplus
}
#endif