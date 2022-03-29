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

#include "sets.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * normalize[x] returns the value of <x>_S (the greatest element of S
 * less than or equal to x):
 *
 * <x>_S = max {s | s \in S, s <= x}.
 */
extern int *normalize;

/******************************************************************
 ******************************************************************/

/**
 * Insert an element into the specified set, in the case it does not
 * belong to the set yet.
 *
 * Parameters:
 * element - The element to be inseted.
 * set - The set where the element will be inserted.
 *
 * Return:
 * The size of the set after this insertion.
 */
int
insert (int element, Set set)
{

  if (set.size == 0 || set.points[set.size - 1] < element)
    {

      set.points[set.size] = element;
      return set.size + 1;
    }
  else
    {
      return set.size;
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Create a new set with size 0 and maximum size "maxSize".
 *
 * Parameter:
 * maxSize - The maximum size of the set.
 *
 * Return:
 * set - The set created.
 */
Set
newSet (int maxSize)
{
  Set set;
  set.size = 0;
  set.points = new int[maxSize];
  return set;
}

/******************************************************************
 ******************************************************************/

/**
 * This function is used by qsort() function to compare two integers.
 */
int
compare (const void *n1, const void *n2)
{
  int *number1 = (int *)n1;
  int *number2 = (int *)n2;

  if (*number1 == *number2)
    {
      return 0;
    }
  else
    {
      return (*number1 > *number2);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Construct the raster points sets X' and Y', which are defined as
 *
 * X' = {<L - x>_X | x \in X} U {0}
 * Y' = {<W - y>_Y | y \in Y} U {0}
 *
 * where X and Y are the integer conic combinations sets for L and W,
 * respectively, and <s'>_S = max {s | s \in S, s <= s'}.
 *
 * Parameters:
 * L                 - Length of the rectangle.
 * W                 - Width of the rectangle.
 * rasterPointsX     - Pointer to the raster points set X'.
 * rasterPointsY     - Pointer to the raster points set Y'.
 * conicCombinations - Set of integer conic combinations of l and w.
 *
 * Remark: it suposes that the integer conic combinations set is
 * sorted.
 */
void
constructRasterPoints (int L, int W, Set *rasterPointsX, Set *rasterPointsY,
                       Set conicCombinations)
{

  int x, i, xSize, ySize;

  /* Maximum raster points X size. */
  xSize = conicCombinations.size;
  for (i = xSize - 1; i >= 0 && conicCombinations.points[i] > L; i--)
    xSize--;

  /* Maximum raster points Y size. */
  ySize = conicCombinations.size;
  for (i = ySize - 1; i >= 0 && conicCombinations.points[i] > W; i--)
    ySize--;

  *rasterPointsX = newSet (xSize + 2);
  *rasterPointsY = newSet (ySize + 2);

  /* Construct the raster points for L. */
  for (i = xSize - 1; i >= 0; i--)
    {
      x = normalize[L - conicCombinations.points[i]];
      (*rasterPointsX).size = insert (x, (*rasterPointsX));
    }

  /* Construct the raster points for W. */
  for (i = ySize - 1; i >= 0; i--)
    {
      x = normalize[W - conicCombinations.points[i]];
      (*rasterPointsY).size = insert (x, (*rasterPointsY));
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Construct the set X of integer conic combinations of l and w.
 *
 * X = {x | x = rl + sw, x <= L, r,s >= 0 integers}
 *
 * Parameters:
 * L - Length of the rectangle.
 * l - Length of the boxes to be packed.
 * w - Width of the boxes to be packed.
 * X - Pointer to the set X.
 */
void
constructConicCombinations (int L, int l, int w, Set *X)
{

  int *inX = (int *)calloc (L + 2, sizeof (int));
  int *c = (int *)calloc (L + 2, sizeof (int));

  (*X).points[0] = 0;
  (*X).size = 1;
  inX[0] = 1;

  for (int i = 0; i <= L; i++)
    {
      c[i] = 0;
    }

  for (int i = l; i <= L; i++)
    {
      if (c[i] < c[i - l] + l)
        {
          c[i] = c[i - l] + l;
        }
    }

  for (int i = w; i <= L; i++)
    {
      if (c[i] < c[i - w] + w)
        {
          c[i] = c[i - w] + w;
        }
    }

  for (int i = 1; i <= L; i++)
    {
      if (c[i] == i && inX[i] == 0)
        {
          (*X).points[(*X).size] = i;
          (*X).size++;
          inX[i] = 1;
        }
    }

  if ((*X).points[(*X).size - 1] != L)
    {
      /* Insert L int the set X. */
      (*X).points[(*X).size] = L;
      (*X).size++;
    }

  free (inX);
  free (c);
}
