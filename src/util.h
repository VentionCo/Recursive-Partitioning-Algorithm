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

#ifndef UTIL_H_
#define UTIL_H_

#define HOMOGENEOUS 0
#define B1 1
#define B2 2
#define B3 3
#define B4 4
#define B5 5
#define B6 6
#define B7 7
#define B8 8
#define B9 9

#define MEM_TYPE_1 1 /* 1-dimensional array. */
#define MEM_TYPE_2 2 /* 2-dimensional array. */
#define MEM_TYPE_3 3 /* 3-dimensional array. */
#define MEM_TYPE_4 4 /* 4-dimensional array. */

#define HORIZONTAL 0
#define VERTICAL 1

#define HORIZONTAL_CUT 0
#define VERTICAL_CUT 1

struct CutPoint
{
  int x1, x2, y1, y2, homogeneous;
};

const unsigned int ptoDiv1 = 2047;
const unsigned int ptoDiv2 = 4192256;
const unsigned int ptoDiv3 = 4290772992UL;

const int nRet = 134217727;
const int solucao = 2013265920;
const int descSol = 27;

const short descPtoDiv2 = 11;
const short descPtoDiv3 = 22;

/******************************************************************
 ******************************************************************/

/*
 * Normalize the L-piece (q0,q1,q2,q3) considering the symmetries
 * below, defined in [3].
 *
 * [3] R. Morabito and S. Morales. A simple and effective recursive
 *     procedure for the manufacturer's pallet loading
 *     problem. Journal of the Operational Research Society, volume
 *     49, number 8, pages 819-828, 1998.
 *
 * Symmetry considerations for (i,j,i',j'):
 *
 * (1) i >= i' and j >= j', from the definition of standardly positioned
 * L's;
 *
 * (2) i >= j, otherwise we could use (j, i, j', i');
 *
 * (3) i = j implies i' >= j', otherwise we again could use (j,i,j',i');
 *
 * (4) i = i' if and only if j = j'. This equivalence follows to avoid
 * degenerated L's which are not explicit rectangles: in terms of
 * occupancy, (i,j,i,j') with j' < j can be replaced by (i,j,i,j) and
 * (i,j,i',j) with i' < i can also be replaced by (i,j,i,j).
 *
 * The normalization (i,j,i',j')^N of a quadruple (i,j,i',j') is defined
 * as:
 *
 * - if (i,j,i',j') satisfies (1)-(4), then (i,j,i',j')^N = (i,j,i',j');
 *
 * - if 0 = i' < i and 0 < j' < j, then (i,j,i',j')^N = (i,j',i,j');
 *
 * - if 0 < i' < i and 0 = j' < j, then (i,j,i',j')^N = (i',j,i',j);
 *
 * - if 0 < i' = i and 0 < j' < j, then (i,j,i',j')^N = (i,j,i,j);
 *
 * - if 0 < i' < i and 0 < j' = j, then (i,j,i',j')^N = (i,j,i,j);
 *
 * - if 0 < i' < i, 0 < j' < j and i < j, then (i,j,i',j')^N = (j,i,j',i');
 *
 * - if 0 < i' < i, 0 < j' < j, i = j and i' < j', then (i,j,i',j')^N =
 *   (j,i,j',i');
 *
 *
 * Parameters:
 * q - The L-piece to be normalized.
 */
void normalizePiece (int *q);

/******************************************************************
 ******************************************************************/

/**
 * Return the index associated to the L-shaped piece (q0, q1, q2, q3).
 */
int getKey (int q0, int q1, int q2, int q3, int type);

/******************************************************************
 ******************************************************************/

/**
 * Return the index associated to the L-shaped piece (q0, q1, q2, q3).
 */
int LIndex (int q0, int q1, int q2, int q3, int memory_type);

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-shaped piece in two new L-shaped pieces, according to
 * the subdivision B1, and put them in the standard position.
 *
 *                  (X,Y)
 * +------------+     o
 * |            |
 * |            |(x,y)                     (x,Y-y')                     (X,y)
 * |      +-----o-----+         +------+     o         +-----------+      o
 * |  L1  |           |         |      |               |           |
 * |      |     L2    |   -->   |      |(x',Y-y)       |           |(X-x',y')
 * +------o           |         |  L1  o-----+         |     L2    o------+
 * |   (x',y')        |         |            |         |                  |
 * |                  |         |            |         |                  |
 * +------------------+         +------------+         +------------------+
 *
 * Parameters:
 * i  - Array of two elements such that i[0] = x' and i[1] = y'.
 *
 * q  - The L-shaped piece to be divided. q = {X, Y, x, y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB1 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-shaped piece in two new L-shaped pieces, according to
 * the subdivision B2, and put them in the standard position.
 *
 *                  (X,Y)
 * +------------+     o
 * |            |
 * |   (x',y')  |                          (x,Y-y)                     (X,y')
 * +------o     |               +-----+      o        +------+           o
 * |      | L1  |               |     |               |      |
 * |      |     |(x,y)    -->   |     |(x-x',Y-y')    |      |(x',y)
 * |      +-----o-----+         |     o------+        |      o-----------+
 * |  L2              |         |  L1        |        |  L2              |
 * |                  |         |            |        |                  |
 * +------------------+         +------------+        +------------------+
 *
 * Parameters:
 * i  - Array of two elements such that i[0] = x' and i[1] = y'.
 *
 * q  - The L-shaped piece to be divided. q = {X, Y, x, y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB2 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-shaped piece in two new L-shaped pieces, according to
 * the subdivision B3, and put them in the standard position.
 *
 *                  (X,Y)                         (X,Y)
 * +------+-----+     o         +------+           o
 * |      |     |               |      |
 * |      |     |(x,y)          |      |                         (X-x',Y-y')
 * |      | L2  o-----+         |      |                  +-----+     o
 * |      |           |         |      |                  |     |
 * |  L1  |           |   -->   |  L1  |(x',y')           |     |(x-x',y-y')
 * |      o-----------+         |      o-----------+      | L2  o-----+
 * |   (x',y')        |         |                  |      |           |
 * |                  |         |                  |      |           |
 * +------------------+         +------------------+      +-----------+
 *
 * Parameters:
 * i  - Array of two elements such that i[0] = x' and i[1] = y'.
 *
 * q  - The L-shaped piece to be divided. q = {X, Y, x, y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB3 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-shaped piece in two new L-shaped pieces, according to
 * the subdivision B4, and put them in the standard position.
 *
 *                  (X,Y)                 (x',Y)
 * +------+           o         +------+     o
 * |      |                     |      |
 * |      |(x,y)                |      |                     (X-x,y)
 * |      o-----------+         |      |             +-----+     o
 * |  L1  |           |         |  L1  |             |     |
 * |      |  (x',y')  |   -->   |      |(x,y')       |     |(X-x',y-y')
 * |      +-----o     |         |      o-----+       |     o-----+
 * |            | L2  |         |            |       |  L2       |
 * |            |     |         |            |       |           |
 * +------------+-----+         +------------+       +-----------+
 *
 * Parameters:
 * i  - Array of two elements such that i[0] = x' and i[1] = y'.
 *
 * q  - The L-shaped piece to be divided. q = {X, Y, x, y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB4 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-shaped piece in two new L-shaped pieces, according to
 * the subdivision B5, and put them in the standard position.
 *
 *                  (X,Y)                  (x,Y)
 * +------------+     o         +------+     o
 * |            |               |      |
 * |     L1     |(x,y)          |      |                    (X-x',y)
 * |            o-----+         |      |(x',Y-y')    +-----+     o
 * |   (x',y')  |     |         |      o-----+       |     |
 * |      o-----+     |   -->   |            |       |     o-----+
 * |      |           |         |     L1     |       | (X-x,y')  |
 * |      |     L2    |         |            |       |           |
 * |      |           |         |            |       |    L2     |
 * +------+-----------+         +------------+       +-----------+
 *
 * Parameters:
 * i  - Array of two elements such that i[0] = x' and i[1] = y'.
 *
 * q  - The L-shaped piece to be divided. q = {X, Y, x, y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB5 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the rectangle in two L-shaped pieces, according to the
 * subdivision B6, and put them in the standard position.
 *
 *                      (X,Y)                 (x'',Y)              (X-x',Y)
 * +-------------+--------o         +------+      o     +--------+      o
 * |             |        |         |      |            |        |
 * |   (x',y')   |   L2   |         |      |            |        |(X-x'',y')
 * |      o------o        |   -->   |      |(x',Y-y')   |        o------+
 * |      |  (x'',y')     |         |  L1  o------+     |   L2          |
 * |  L1  |               |         |             |     |               |
 * |      |               |         |             |     |               |
 * +------+---------------+         +-------------+     +---------------+
 *
 * Parameters:
 * i  - Array of three elements such that i[0] = x', i[1] = y' and i[2] = x''.
 *
 * q  - The rectangle to be divided. q = {X, Y, X, Y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */

void standardPositionB6 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the rectangle in two L-shaped pieces, according to the
 * subdivision B7, and put them in the standard position.
 *
 *             (X,Y)
 * +-------------o
 * |             |
 * |   (x',y'')  |                                           (X,y'')
 * |      o------+                     (X,Y-y')    +------+      o
 * |      |      |         +------+      o         |      |
 * |  L1  |  L2  |         |      |                |      |
 * |      |      |   -->   |      |                |      |(X-x',y')
 * +------o      |         |      |(x',Y-y'')      |      o------+
 * |   (x',y')   |         |  L1  o------+         |             |
 * |             |         |             |         |      L2     |
 * |             |         |             |         |             |
 * +-------------+         +-------------+         +-------------+
 *
 * Parameters:
 * i  - Array of three elements such that i[0] = x', i[1] = y' and i[2] = y''.
 *
 * q  - The rectangle to be divided. q = {X, Y, X, Y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB7 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-shaped piece in two new L-shaped pieces, according to
 * the subdivision B8, and put them in the standard position.
 *
 *                  (X,Y)                  (x,Y)
 * +------------+     o         +------+     o
 * |            |               |      |
 * |   (x',y')  |               |      |                     (X-x',y')
 * |      o-----+               |      |              +-----+     o
 * |      |     |               |  L1  |              |     |
 * |  L1  |     |(x,y)    -->   |      |(x',Y-y')     |     |(x-x',y)
 * |      |     o-----+         |      o-----+        |     o-----+
 * |      |  L2       |         |            |        |  L2       |
 * |      |           |         |            |        |           |
 * +------+-----------+         +------------+        +-----------+
 *
 * Parameters:
 * i  - Array of two elements such that i[0] = x' and i[1] = y'.
 *
 * q  - The L-shaped piece to be divided. q = {X, Y, x, y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB8 (int *i, int *q, int *q1, int *q2);

/******************************************************************
 ******************************************************************/

/**
 * Divide the L-shaped piece in two new L-shaped pieces, according to
 * the subdivision B9, and put them in the standard position.
 *
 *                  (X,Y)
 * +---------+        o
 * |         |
 * |         |(x,y)                                                      (X,y)
 * |   L1    o---+----+                    (x',Y-y')    +----+             o
 * |             |    |         +---------+   o         |    |
 * |             |    |         |         |             |    |(X-x',y')
 * +-------------o    |   -->   |         |(x,y-y')     |    o-------------+
 * |          (x',y') |         |   L1    o---+         |                  |
 * |     L2           |         |             |         |        L2        |
 * |                  |         |             |         |                  |
 * +------------------+         +-------------+         +------------------+
 *
 * Parameters:
 * i  - Array of two elements such that i[0] = x' and i[1] = y'.
 *
 * q  - The L-shaped piece to be divided. q = {X, Y, x, y}.
 *
 * q1 - Array to store L1.
 *
 * q2 - Array to store L2.
 */
void standardPositionB9 (int *i, int *q, int *q1, int *q2);

#endif
