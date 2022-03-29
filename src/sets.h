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

#ifndef SETS_H_
#define SETS_H_

struct Set
{
  int size;
  int *points;
};

/**
 * Create a new set with size 0 and maximum size "maxSize".
 *
 * Parameter:
 * maxSize - The maximum size of the set.
 *
 * Return:
 * set - The set created.
 */
Set newSet (int maxSize);

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
void constructRasterPoints (int L, int W, Set *rasterPointsX,
                            Set *rasterPointsY, Set conicCombinations);

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
void constructConicCombinations (int L, int l, int w, Set *X);

#endif
