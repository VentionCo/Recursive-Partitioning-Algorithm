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

#ifndef BD_H_
#define BD_H_

/**
 * Guillotine and first order non-guillotine cuts recursive procedure.
 *
 * Parameters:
 * L     - Length of the pallet.
 * W     - Width of the pallet.
 * l     - Length of the boxes.
 * w     - Width of the boxes.
 * N_max - Maximum search depth.
 *
 * Return:
 * - the number of (l,w)-boxes packed into (L,W) pallet.
 */
int solve_BD (int L, int W, int l, int w, int N_max);

#endif
