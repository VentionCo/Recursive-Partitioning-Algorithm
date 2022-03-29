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

#include "draw_bd.h"
#include "util.h"

#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>

extern int l, w;

extern const int **lowerBound, *indexX, *indexY;

extern const int *divisionPoint;
extern const int *solution;

extern std::map<int, int> *solutionMap;
extern std::map<int, int> *divisionPointMap;

extern const int *normalize;
extern const int memory_type;

int ret;
int **ptoRet;

void drawR (int L, int *q);

/******************************************************************
 ******************************************************************/

inline int
LIndex (int q0, int q1, int q2, int q3)
{
  return LIndex (q0, q1, q2, q3, memory_type);
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
 * The lower bound.
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
 * Determine how to cut the L-piece.
 *
 * Parameter:
 * q - The L-piece.
 */
short
LCut (int *q)
{
  /* Divide the L-piece in two rectangles and calculate their lower
   * bounds to compose the lower bound of the L-piece. */
  int a = R_LowerBound (q[2], q[1]) + R_LowerBound (q[0] - q[2], q[3]);
  int b = R_LowerBound (q[2], q[1] - q[3]) + R_LowerBound (q[0], q[3]);

  return (a > b) ? VERTICAL_CUT : HORIZONTAL_CUT;
}

/******************************************************************
 ******************************************************************/

/**
 * Fix the coordinates of the rectangle with the specified identifier.
 *
 * Paremeter:
 * id - Identifier of the rectangle.
 */
void
fixCoordinates (int id)
{
  if (ptoRet[id][0] > ptoRet[id][2])
    {
      std::swap (ptoRet[id][0], ptoRet[id][2]);
    }
  if (ptoRet[id][1] > ptoRet[id][3])
    {
      std::swap (ptoRet[id][1], ptoRet[id][3]);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Normalize degenerated L-pieces.
 *
 * Parameter:
 * q - The degenerated L-piece to be normalized.
 */
void
normalizeDegeneratedL (int *q)
{
  if (q[2] == 0)
    {
      q[2] = q[0];
      q[1] = q[3];
    }
  else if (q[3] == 0)
    {
      q[3] = q[1];
      q[0] = q[2];
    }
  else if (q[2] == q[0] || q[3] == q[1])
    {
      q[2] = q[0];
      q[3] = q[1];
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Shift the rectangle in the x-axis.
 *
 * Paremeters:
 * id     - Identifier of the rectangle.
 *
 * deltaX - Amount to be shifted.
 */
void
shiftX (int id, int deltaX)
{
  ptoRet[id][0] += deltaX;
  ptoRet[id][2] += deltaX;
}

/******************************************************************
 ******************************************************************/

/**
 * Shift the rectangle in the y-axis.
 *
 * Paremeters:
 * id     - Identifier of the rectangle.
 *
 * deltaY - Amount to be shifted.
 */
void
shiftY (int id, int deltaY)
{
  ptoRet[id][1] += deltaY;
  ptoRet[id][3] += deltaY;
}

/******************************************************************
 ******************************************************************/

void
P1 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      ptoRet[i][1] = q[1] - ptoRet[i][1];
      ptoRet[i][3] = q[1] - ptoRet[i][3];

      fixCoordinates (i);

      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
P2 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      ptoRet[i][0] = q[0] - ptoRet[i][0];
      ptoRet[i][2] = q[0] - ptoRet[i][2];

      fixCoordinates (i);

      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
P3 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      ptoRet[i][0] = q[0] - ptoRet[i][0];
      ptoRet[i][2] = q[0] - ptoRet[i][2];

      ptoRet[i][1] = q[1] - ptoRet[i][1];
      ptoRet[i][3] = q[1] - ptoRet[i][3];

      fixCoordinates (i);

      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
P4 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
P5 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      int tmp1 = ptoRet[i][1];
      int tmp2 = ptoRet[i][3];

      ptoRet[i][1] = q[0] - ptoRet[i][0];
      ptoRet[i][3] = q[0] - ptoRet[i][2];

      ptoRet[i][0] = tmp1;
      ptoRet[i][2] = tmp2;

      fixCoordinates (i);

      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
P6 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      int tmp1 = ptoRet[i][0];
      int tmp2 = ptoRet[i][2];

      ptoRet[i][0] = q[1] - ptoRet[i][1];
      ptoRet[i][2] = q[1] - ptoRet[i][3];

      ptoRet[i][1] = tmp1;
      ptoRet[i][3] = tmp2;

      fixCoordinates (i);

      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
P7 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      int tmp1 = q[0] - ptoRet[i][0];
      int tmp2 = q[0] - ptoRet[i][2];

      ptoRet[i][0] = q[1] - ptoRet[i][1];
      ptoRet[i][2] = q[1] - ptoRet[i][3];

      ptoRet[i][1] = tmp1;
      ptoRet[i][3] = tmp2;

      fixCoordinates (i);

      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
P8 (int start, int end, int L, int *q, int deltaX, int deltaY)
{

  for (int i = start; i < end; i++)
    {
      int tmp1 = ptoRet[i][0];
      int tmp2 = ptoRet[i][2];

      ptoRet[i][0] = ptoRet[i][1];
      ptoRet[i][2] = ptoRet[i][3];

      ptoRet[i][1] = tmp1;
      ptoRet[i][3] = tmp2;

      shiftX (i, deltaX);
      shiftY (i, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B1 subdivision.
 */
void
drawB1 (int L, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[2];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L] & ptoDiv1;
      div[1] = (divisionPoint[L] & ptoDiv2) >> descPtoDiv2;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L][h] & ptoDiv1;
      div[1] = (divisionPointMap[L][h] & ptoDiv2) >> descPtoDiv2;
    }

  standardPositionB1 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = div[1];
  if (div[0] == 0)
    {
      deltaY = q[3];
    }

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L1, q1, deltaX, deltaY);
      else
        P8 (start, end, L1, q1, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P1 (start, end, L1, q1, deltaX, deltaY);
      else
        P5 (start, end, L1, q1, deltaX, deltaY);
    }

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = 0;
  deltaY = 0;
  if (div[1] == 0)
    deltaX = div[0];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L2, q2, deltaX, deltaY);
      else
        P8 (start, end, L2, q2, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P2 (start, end, L2, q2, deltaX, deltaY);
      else
        P6 (start, end, L2, q2, deltaX, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B2 subdivision.
 */
void
drawB2 (int L, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[3];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L] & ptoDiv1;
      div[1] = (divisionPoint[L] & ptoDiv2) >> descPtoDiv2;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L][h] & ptoDiv1;
      div[1] = (divisionPointMap[L][h] & ptoDiv2) >> descPtoDiv2;
    }

  standardPositionB2 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = q[3];
  if (div[1] == q[1])
    deltaX = div[0];
  else if (tmp[0] == tmp[2])
    deltaY = div[1];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L1, q1, deltaX, deltaY);
      else
        P8 (start, end, L1, q1, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P3 (start, end, L1, q1, deltaX, deltaY);
      else
        P7 (start, end, L1, q1, deltaX, deltaY);
    }

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = 0;
  deltaY = 0;

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (width >= height)
    P4 (start, end, L2, q2, deltaX, deltaY);
  else
    P8 (start, end, L2, q2, deltaX, deltaY);
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B3 subdivision.
 */
void
drawB3 (int L, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[2];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L] & ptoDiv1;
      div[1] = (divisionPoint[L] & ptoDiv2) >> descPtoDiv2;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L][h] & ptoDiv1;
      div[1] = (divisionPointMap[L][h] & ptoDiv2) >> descPtoDiv2;
    }

  standardPositionB3 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = 0;

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (width >= height)
    P4 (start, end, L1, q1, deltaX, deltaY);
  else
    P8 (start, end, L1, q1, deltaX, deltaY);

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = div[0];
  deltaY = div[1];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (width >= height)
    P4 (start, end, L2, q2, deltaX, deltaY);
  else
    P8 (start, end, L2, q2, deltaX, deltaY);
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B4 subdivision.
 */
void
drawB4 (int L, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[2];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L] & ptoDiv1;
      div[1] = (divisionPoint[L] & ptoDiv2) >> descPtoDiv2;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L][h] & ptoDiv1;
      div[1] = (divisionPointMap[L][h] & ptoDiv2) >> descPtoDiv2;
    }

  standardPositionB4 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = 0;

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (width >= height)
    P4 (start, end, L1, q1, deltaX, deltaY);
  else
    P8 (start, end, L1, q1, deltaX, deltaY);

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = q[3];
  deltaY = 0;
  if (div[0] == q[0])
    deltaY = div[1];
  else if (tmp[0] == tmp[2])
    deltaX = div[0];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L2, q2, deltaX, deltaY);
      else
        P8 (start, end, L2, q2, deltaX, deltaY);
    }
  else
    {
      if (width > height)
        P3 (start, end, L2, q2, deltaX, deltaY);
      else
        P7 (start, end, L2, q2, deltaX, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B5 subdivision.
 */
void
drawB5 (int L, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[2];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L] & ptoDiv1;
      div[1] = (divisionPoint[L] & ptoDiv2) >> descPtoDiv2;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L][h] & ptoDiv1;
      div[1] = (divisionPointMap[L][h] & ptoDiv2) >> descPtoDiv2;
    }

  standardPositionB5 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = 0;
  if (div[0] == 0)
    deltaY = div[1];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L1, q1, deltaX, deltaY);
      else
        P8 (start, end, L1, q1, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P1 (start, end, L1, q1, deltaX, deltaY);
      else
        P5 (start, end, L1, q1, deltaX, deltaY);
    }

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = div[0];
  deltaY = 0;
  if (div[1] == 0)
    deltaX = q[2];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L2, q2, deltaX, deltaY);
      else
        P8 (start, end, L2, q2, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P2 (start, end, L2, q2, deltaX, deltaY);
      else
        P6 (start, end, L2, q2, deltaX, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B6 subdivision.
 */
void
drawB6 (int L, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[3];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L] & ptoDiv1;
      div[1] = (divisionPoint[L] & ptoDiv2) >> descPtoDiv2;
      div[2] = (divisionPoint[L] & ptoDiv3) >> descPtoDiv3;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L][h] & ptoDiv1;
      div[1] = (divisionPointMap[L][h] & ptoDiv2) >> descPtoDiv2;
      div[2] = (divisionPointMap[L][h] & ptoDiv3) >> descPtoDiv3;
    }

  standardPositionB6 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = 0;
  if (div[0] == 0)
    deltaY = div[1];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L1, q1, deltaX, deltaY);
      else
        P8 (start, end, L1, q1, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P1 (start, end, L1, q1, deltaX, deltaY);
      else
        P5 (start, end, L1, q1, deltaX, deltaY);
    }

  /*Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = div[0];
  deltaY = 0;
  if (div[1] == 0)
    deltaX = div[2];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L2, q2, deltaX, deltaY);
      else
        P8 (start, end, L2, q2, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P2 (start, end, L2, q2, deltaX, deltaY);
      else
        P6 (start, end, L2, q2, deltaX, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B7 subdivision.
 */
void
drawB7 (int L, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[3];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L] & ptoDiv1;
      div[1] = (divisionPoint[L] & ptoDiv2) >> descPtoDiv2;
      div[2] = (divisionPoint[L] & ptoDiv3) >> descPtoDiv3;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L][h] & ptoDiv1;
      div[1] = (divisionPointMap[L][h] & ptoDiv2) >> descPtoDiv2;
      div[2] = (divisionPointMap[L][h] & ptoDiv3) >> descPtoDiv3;
    }

  standardPositionB7 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = div[1];
  if (div[0] == 0)
    deltaY = div[2];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L1, q1, deltaX, deltaY);
      else
        P8 (start, end, L1, q1, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P1 (start, end, L1, q1, deltaX, deltaY);
      else
        P5 (start, end, L1, q1, deltaX, deltaY);
    }

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = 0;
  deltaY = 0;
  if (div[1] == 0)
    deltaX = div[0];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L2, q2, deltaX, deltaY);
      else
        P8 (start, end, L2, q2, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P2 (start, end, L2, q2, deltaX, deltaY);
      else
        P6 (start, end, L2, q2, deltaX, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B8 subdivision.
 */
void
drawB8 (int L_index, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[2];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L_index] & ptoDiv1;
      div[1] = (divisionPoint[L_index] & ptoDiv2) >> descPtoDiv2;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L_index][h] & ptoDiv1;
      div[1] = (divisionPointMap[L_index][h] & ptoDiv2) >> descPtoDiv2;
    }

  standardPositionB8 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = 0;
  if (div[0] == 0)
    deltaY = div[1];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L1, q1, deltaX, deltaY);
      else
        P8 (start, end, L1, q1, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P1 (start, end, L1, q1, deltaX, deltaY);
      else
        P5 (start, end, L1, q1, deltaX, deltaY);
    }

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = div[0];
  deltaY = 0;
  if (div[1] == 0)
    deltaX = div[0];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (width >= height)
    P4 (start, end, L2, q2, deltaX, deltaY);
  else
    P8 (start, end, L2, q2, deltaX, deltaY);
}

/******************************************************************
 ******************************************************************/

/**
 * Draw the boxes according to the B9 subdivision.
 */
void
drawB9 (int L_index, int *q)
{
  int L1, L2;
  int q1[4], q2[4], tmp[4];
  int width, height;
  int deltaX, deltaY;
  int start, end;
  int div[2];

  if (memory_type == MEM_TYPE_4)
    {
      div[0] = divisionPoint[L_index] & ptoDiv1;
      div[1] = (divisionPoint[L_index] & ptoDiv2) >> descPtoDiv2;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      div[0] = divisionPointMap[L_index][h] & ptoDiv1;
      div[1] = (divisionPointMap[L_index][h] & ptoDiv2) >> descPtoDiv2;
    }

  standardPositionB9 (div, q, q1, q2);

  /* Draw L1. */
  tmp[0] = q1[0];
  tmp[1] = q1[1];
  tmp[2] = q1[2];
  tmp[3] = q1[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q1);
  L1 = LIndex (q1[0], q1[1], q1[2], q1[3]);

  start = ret;
  drawR (L1, q1);
  end = ret;

  deltaX = 0;
  deltaY = div[1];
  if (div[0] == 0)
    deltaY = q[3];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (width >= height)
    P4 (start, end, L1, q1, deltaX, deltaY);
  else
    P8 (start, end, L1, q1, deltaX, deltaY);

  /* Draw L2. */
  tmp[0] = q2[0];
  tmp[1] = q2[1];
  tmp[2] = q2[2];
  tmp[3] = q2[3];
  normalizeDegeneratedL (tmp);

  normalizePiece (q2);
  L2 = LIndex (q2[0], q2[1], q2[2], q2[3]);

  start = ret;
  drawR (L2, q2);
  end = ret;

  deltaX = 0;
  deltaY = 0;
  if (div[1] == 0)
    deltaX = div[0];

  if (tmp[0] != tmp[1])
    {
      width = tmp[0];
      height = tmp[1];
    }
  else
    {
      width = tmp[2];
      height = tmp[3];
    }

  if (tmp[0] == tmp[2])
    {
      if (width >= height)
        P4 (start, end, L2, q2, deltaX, deltaY);
      else
        P8 (start, end, L2, q2, deltaX, deltaY);
    }
  else
    {
      if (width >= height)
        P2 (start, end, L2, q2, deltaX, deltaY);
      else
        P6 (start, end, L2, q2, deltaX, deltaY);
    }
}

/******************************************************************
 ******************************************************************/

void
drawR (int L, int *q)
{
  int i;
  int start, end;
  int divisionType;

  if (memory_type == MEM_TYPE_4)
    {
      divisionType = (solution[L] & solucao) >> descSol;
    }
  else
    {
      int h = getKey (q[0], q[1], q[2], q[3], memory_type);
      divisionType = (solutionMap[L][h] & solucao) >> descSol;
    }

  switch (divisionType)
    {

    case HOMOGENEOUS:

      /* Non-degenerated L. */
      if (q[0] != q[2])
        {

          short cut = LCut (q);

          if (cut == VERTICAL_CUT)
            {
              ret = drawBD (q[2], q[1], ret);
              start = ret;
              ret = drawBD (normalize[q[0] - q[2]], q[3], ret);
              end = ret;
              for (i = start; i < end; i++)
                {
                  shiftX (i, q[2]);
                }
            }
          else
            {
              start = ret;
              ret = drawBD (q[2], normalize[q[1] - q[3]], ret);
              end = ret;
              for (i = start; i < end; i++)
                {
                  shiftY (i, q[3]);
                }
              ret = drawBD (q[0], q[3], ret);
            }
        }
      /* Degenerated L (rectangle). */
      else
        {
          ret = drawBD (q[0], q[1], ret);
        }
      break;

    case B1:
      drawB1 (L, q);
      break;
    case B2:
      drawB2 (L, q);
      break;
    case B3:
      drawB3 (L, q);
      break;
    case B4:
      drawB4 (L, q);
      break;
    case B5:
      drawB5 (L, q);
      break;
    case B6:
      drawB6 (L, q);
      break;
    case B7:
      drawB7 (L, q);
      break;
    case B8:
      drawB8 (L, q);
      break;
    case B9:
      drawB9 (L, q);
      break;
    default:
      break;
    }
}

/******************************************************************
 ******************************************************************/

std::string
MakeJsonString (int Lo, int Wo, int L, int *q, int n, int l, int w, bool swap)
{
  std::string str = "[";;
  float x, y;
  float xl, yl, xh, yh;
  bool rotated = true;
  int sym = 0 == l - w;
  int cmp = swap ? l : w;
  if (sym)
    rotated = false;

  printf ("[");
  for (int i = 0; i < n; i++) {
    if (!sym)
      {
        if (sym || 0 == ptoRet[i][3] - ptoRet[i][1] - cmp)
          rotated = false;
        else
          rotated = true;
      }
    xh = (float)ptoRet[i][0];
    yh = (float)ptoRet[i][1];
    xl = (float)ptoRet[i][2];
    yl = (float)ptoRet[i][3];
    x = (xl - xh) / 2.0 + xh;
    y = (yl - yh) / 2.0 + yh;

    str += "{\"x\": " +  std::to_string(swap ? x : y) + ", \"y\": " + std::to_string(swap ? y : x) + ", \"rotated\": " + (rotated ? "true" : "false") + "}";
    if (i + 1 != n) {
      str += ",\n";
    }
  }

  str += "]";
  return str;
}

/******************************************************************
 ******************************************************************/

std::string
draw (int Lo, int Wo, int L, int *q, int n, bool solvedWithL, int l, int w,
      bool swap)
{
  ptoRet = (int **)malloc ((n) * sizeof (int *));
  if (ptoRet == NULL)
    {
      printf ("Error allocating memory.\n");
      exit (0);
    }

  for (int i = 0; i < n; i++)
    {
      ptoRet[i] = (int *)malloc (4 * sizeof (int));
      if (ptoRet[i] == NULL)
        {
          printf ("Error allocating memory.\n");
          exit (0);
        }
    }

  ret = 0;
  if (solvedWithL)
    {
      drawR (L, q);
    }
  else
    {
      ret = drawBD (q[0], q[1], ret);
    }

  std::string result = MakeJsonString (Lo, Wo, L, q, n, l, w, swap);

  // for (int i = 0; i < n; i++)
  //   free (ptoRet[i]);
  free (ptoRet);

  return result;
}
