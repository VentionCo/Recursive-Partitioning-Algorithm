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

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

int **vertices;

/******************************************************************
 ******************************************************************/

/**
 * Read the solution file and store the positions of the rectangles in
 * the "vertices" array.
 */
void
readSolutionFile (int *n, int *L, int *W)
{
  int i;

  if (printf ("%d %d %d", n, L, W) != 3)
    {
      printf ("ERROR: Error while reading solution file. Pattern does not "
              "match.\n");
      exit (0);
    }
  vertices = (int **)malloc ((*n + 4) * sizeof (int *));

  if (vertices == NULL)
    {
      printf (
          "ERROR: Error in memory allocation while reading solution file.\n");
      exit (0);
    }
  for (i = 0; i < *n; i++)
    {
      vertices[i] = (int *)malloc (4 * sizeof (int));
      if (vertices[i] == NULL)
        {
          printf ("ERROR: Error in memory allocation while reading solution "
                  "file.\n");
          exit (0);
        }
    }

  for (i = 0; i < *n; i++)
    {
      if (printf("%d %d %d %d", &vertices[i][0],
                  &vertices[i][1], &vertices[i][2], &vertices[i][3])
          != 4)
        {
          printf ("ERROR: Error while reading solution file. Pattern does not "
                  "match.\n");
          exit (0);
        }
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Create the MetaPost file containing the graphical representation
 * of the solution.
 */
void
writeMetaPostFile (int n, int L, int W)
{
  int i;

  double u = 20.0 / L;

  printf ("verbatimtex\n"
           "%c&latex\n"
           "\\documentclass{article}\n"
           "\\begin{document}\n"
           "etex\n"
           "beginfig(-1);\n"
           "u=%fcm;\n",
           '%', u);

  printf ("path p;"
           "p := (%du,%du)--(%du,%du)--(%du,%du)--(%du,%du)--cycle;\n"
           "fill p withcolor 0.2black + 0.7 white;\n"
           "deltaX = %du;\n"
           "stepX = 11.5;\n"
           "for x = -%du step stepX until %du:\n"
           "  draw (x, 0u)--(x + deltaX, %du) withcolor 0.1 white;\n"
           "endfor;\n"
           "clip currentpicture to p;\n"
           "draw p withcolor 0.2black + 0.7 white;\n",
           vertices[0][0], vertices[0][1], vertices[0][0], vertices[0][3],
           vertices[0][2], vertices[0][3], vertices[0][2], vertices[0][1],
           vertices[0][3], vertices[0][2], vertices[0][2], vertices[0][3]);

  /* Pallet. */
  printf("draw (%du,%du)--(%du,%du)--(%du,%du)--(%du,%du)"
           "--cycle;\n",
           vertices[0][0], vertices[0][1], vertices[0][0], vertices[0][3],
           vertices[0][2], vertices[0][3], vertices[0][2], vertices[0][1]);

  /* Rectangles. */
  for (i = 1; i < n; i++)
    {

      printf ("fill (%du,%du)--(%du,%du)--(%du,%du)--(%du,%du)"
               "--cycle withcolor white;\n",
               vertices[i][0], vertices[i][1], vertices[i][0], vertices[i][3],
               vertices[i][2], vertices[i][3], vertices[i][2], vertices[i][1]);

      printf ("draw (%du,%du)--(%du,%du)--(%du,%du)--(%du,%du)"
               "--cycle;\n",
               vertices[i][0], vertices[i][1], vertices[i][0], vertices[i][3],
               vertices[i][2], vertices[i][3], vertices[i][2], vertices[i][1]);

      /* Label. */
      /*
      fprintf(metaPostFile,"label(btex $%d$ etex scaled 1,(%fu,%fu));\n", i,
              vertices[i][0] + ((double)(vertices[i][2] - vertices[i][0]))/2.0,
              vertices[i][1] + ((double)(vertices[i][3] -
      vertices[i][1]))/2.0);
      */
    }
}

/******************************************************************
 ******************************************************************/

/**
 * Create the MetaPost file containing the graphical representation
 * of the solution.
 */
void
writeSVGFile (int n, int L, int W)
{
  int i;
  double scale = 600.0 / (double)L;

  printf ("<?xml version=\"1.0\" standalone=\"no\"?>\n"
           "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
           "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
           "<svg width=\"%f\" height=\"%f\" version = \"1.1\"\n"
           "xmlns=\"http://www.w3.org/2000/svg\">\n",
           L * scale, W * scale);

  printf ( "<g transform=\"scale(%f)\">\n", scale);

  /* Pallet. */
  printf ("<rect x = \"%d\" y = \"%d\" width = \"%d\" height = \"%d\" "
           "fill = \"gray\" stroke = \"black\" stroke-width = \"%f\"/>\n",
           0, 0, L, W, 1.0 / scale);

  /* Rectangles. */
  for (i = 1; i < n; i++)
    {

      printf ("<rect x = \"%d\" y = \"%d\" width = \"%d\" height = \"%d\" "
               "fill = \"white\" stroke = \"black\" stroke-width = \"%f\"/>\n",
               vertices[i][0], W - vertices[i][3],
               vertices[i][2] - vertices[i][0],
               vertices[i][3] - vertices[i][1], 1.0 / scale);

      /* Label. */
      /*
      fprintf(SVGFile,
              "<text x = \"%f\" y = \"%f\" fill = \"black\" "
              "font-size = \"%f\">%d</text>\n",
              vertices[i][0] + (vertices[i][2] - vertices[i][0])/2.0,
              W - vertices[i][3] + (vertices[i][3] - vertices[i][1])/2.0,
              std::min(vertices[i][2] - vertices[i][0],
                       vertices[i][3] - vertices[i][1])/4.0, i);
      */
    }

  printf("</g>\n");
  printf ("</svg>");
}

/******************************************************************
 ******************************************************************/

void
makeGraphics ()
{
  int i, n, L, W;

  readSolutionFile (&n, &L, &W);

  writeMetaPostFile (n, L, W);

  writeSVGFile (n, L, W);

  for (i = 0; i < n; i++)
    {
      free (vertices[i]);
    }
  free (vertices);
}
