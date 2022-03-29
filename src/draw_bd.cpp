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

#include "util.h"
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

void draw (int L, int W, int dx, int dy);

extern const CutPoint **cutPoints;
extern const int *normalize;
extern const int *indexX, *indexY;
extern int **ptoRet;
extern const int l, w;

int boxesDrawn = 0;

/******************************************************************
 ******************************************************************/

/**
 * Determine the orientation of the boxes (l,w) that maximizes the
 * homogeneous packing in (x,y).
 *
 * Parameters:
 * x - Length of the rectangle.
 *
 * y - Width of the rectangle.
 */
short
boxOrientation (int x, int y)
{
  int a = (x / l) * (y / w);
  int b = (x / w) * (y / l);
  return (a > b) ? HORIZONTAL : VERTICAL;
}

/******************************************************************
 ******************************************************************/

void
drawHomogeneous (int x, int y, int dx, int dy)
{
  int i, j;
  short corte = boxOrientation (x, y);

  if (corte == HORIZONTAL)
    {
      for (i = 0; i + l <= x; i += l)
        {
          for (j = 0; j + w <= y; j += w)
            {
              ptoRet[boxesDrawn][0] = i + dx;
              ptoRet[boxesDrawn][1] = j + dy;
              ptoRet[boxesDrawn][2] = i + l + dx;
              ptoRet[boxesDrawn][3] = j + w + dy;
              boxesDrawn++;
            }
        }
    }

  else
    {
      for (i = 0; i + w <= x; i += w)
        {
          for (j = 0; j + l <= y; j += l)
            {
              ptoRet[boxesDrawn][0] = i + dx;
              ptoRet[boxesDrawn][1] = j + dy;
              ptoRet[boxesDrawn][2] = i + w + dx;
              ptoRet[boxesDrawn][3] = j + l + dy;
              boxesDrawn++;
            }
        }
    }
}

/******************************************************************
 ******************************************************************/

void
getSubproblems (CutPoint cutPoint, int L_[6], int W_[6], int L, int W)
{
  int x1 = cutPoint.x1;
  int x2 = cutPoint.x2;
  int y1 = cutPoint.y1;
  int y2 = cutPoint.y2;

  L_[1] = x1;
  W_[1] = normalize[W - y1];
  L_[2] = normalize[L - x1];
  W_[2] = normalize[W - y2];
  L_[3] = normalize[x2 - x1];
  W_[3] = normalize[y2 - y1];
  L_[4] = x2;
  W_[4] = y1;
  L_[5] = normalize[L - x2];
  W_[5] = y2;
}

/******************************************************************
 ******************************************************************/

void
drawRotation (int L, int W, int dx, int dy)
{
  int L_[6], W_[6];
  int i, iX, iY;

  std::swap (L, W);

  iX = indexX[L];
  iY = indexY[W];

  if (cutPoints[iX][iY].homogeneous)
    {
      std::swap (L, W);
      drawHomogeneous (L, W, dx, dy);
      return;
    }

  getSubproblems (cutPoints[iX][iY], L_, W_, L, W);

  for (i = 1; i <= 5; i++)
    {
      std::swap (L_[i], W_[i]);
      if (L_[i] != 0 && W_[i] != 0
          && !((L_[i] == L && W_[i] == W) || (L_[i] == W && W_[i] == L)))
        {
          switch (i)
            {
            case 1:
              draw (L_[1], W_[1], dx + W_[4], dy + L_[2]);
              break;
            case 2:
              draw (L_[2], W_[2], dx + W_[5], dy);
              break;
            case 3:
              draw (L_[3], W_[3], dx + W_[4], dy + L_[5]);
              break;
            case 4:
              draw (L_[4], W_[4], dx, dy + L_[5]);
              break;
            case 5:
              draw (L_[5], W_[5], dx, dy);
              break;
            }
        }
      std::swap (L_[i], W_[i]);
    }
}

/******************************************************************
 ******************************************************************/

void
drawNormal (int L, int W, int dx, int dy)
{
  int i;
  int L_[6], W_[6];

  int iX = indexX[L];
  int iY = indexY[W];

  if (cutPoints[iX][iY].homogeneous)
    {
      drawHomogeneous (L, W, dx, dy);
      return;
    }

  getSubproblems (cutPoints[iX][iY], L_, W_, L, W);

  for (i = 1; i <= 5; i++)
    {
      if (L_[i] != 0 && W_[i] != 0
          && !((L_[i] == L && W_[i] == W) || (L_[i] == W && W_[i] == L)))
        {
          switch (i)
            {
            case 1:
              draw (L_[1], W_[1], dx, dy + W_[4]);
              break;
            case 2:
              draw (L_[2], W_[2], dx + L_[1], dy + W_[5]);
              break;
            case 3:
              draw (L_[3], W_[3], dx + L_[1], dy + W_[4]);
              break;
            case 4:
              draw (L_[4], W_[4], dx, dy);
              break;
            case 5:
              draw (L_[5], W_[5], dx + L_[4], dy);
              break;
            }
        }
    }
}

/******************************************************************
 ******************************************************************/

void
draw (int L, int W, int dx, int dy)
{
  if (L >= W)
    {
      drawNormal (L, W, dx, dy);
    }
  else
    {
      drawRotation (L, W, dx, dy);
    }
}

int
drawBD (int L, int W, int ret)
{
  boxesDrawn = ret;
  draw (L, W, 0, 0);
  return boxesDrawn;
}
