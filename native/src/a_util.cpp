#include "a_clib.h"
#include "c_vectormath.h"
#include<complex>
#include<cmath>

using namespace std::complex_literals;

// compute the S \cdot (S_i \times S_j)
inline double volume(double S[3], double Si[3], double Sj[3]) {
  double tx = S[0] * (-Si[2] * Sj[1] + Si[1] * Sj[2]);
  double ty = S[1] * (Si[2] * Sj[0] - Si[0] * Sj[2]);
  double tz = S[2] * (-Si[1] * Sj[0] + Si[0] * Sj[1]);
  return tx + ty + tz;
}

double skyrmion_number(double * spin, double * charge, int nx, int ny, int nz,
                       int * ngbs, int n_ngbs) {

  /* Calculation of the "Skyrmion number" Q for a two dimensional discrete
   * spin lattice in the x-y plane (also known
   * as finite spin chirality)
   *
   * The *spin array is the vector field for a two dimensional
   * lattice with dimensions nx * ny
   * (we can take a slice of a bulk from Python and pass it here,
   *  remember to do the ame for the neighbours matrix)
   * The array follows the order:
   *   [Sx0 Sy0 Sz0 Sx1 Sy1 Sz1 ... ]
   *
   * Charge is a scalar field array used to store the spin chirality /
   * skyrmion number density (Q value per lattice site)
   *
   * *ngbs is the array with the neighbours information for every
   * lattice site. The array is like:
   *
   *                                        __ indexes beyond nearest neighbours
   *                                       |
   *      [ 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, ... 1-x, 1+x, 1-y, ...  ]
   *        i=0                               i=n_ngbs       ...
   *
   * where  0-y  is the index of the neighbour of the 0th spin,
   * in the -y direction, for example.
   * The first 6 indexes per every i are the NEAREST neighbours
   *
   *
   * THEORY:
   *
   * Referring to the i-th site of a
   * square lattice, we refer to the x direction for the i-th neighbours
   * and to the y direction for the j-th neighbours, as here:
   *
   *                    X  (j + 1)
   *                    |
   *                (i) |
   *     (i - 1) X------O------X   (i + 1)
   *                    |
   *                    |
   *                    X   (j - 1)
   *
   * Then, we use the following expression:
   *
   * Q =  S_i \dot ( S_{i+1} \times S_{j+1} )
   *      +  S_i \dot ( S_{i-1} \times S_{j-1} )
   *
   * This expression is based on the publication PRL 108, 017601 (2012)
   * where Q is called "finite spin chirality". The idea comes from
   * discrete chiral quantities in Hall effect studies. For example, at
   * the end of page 3 in Rep. Prog. Phys. 78 (2015) 052502, it
   * is argued:
   *     scalar chirality (...) , which measures the volume enclosed
   *     by the three spins of the elementary triangle and, similarly to
   *     (the vector chirlity) is sensitive to the sense of spin's
   *     rotation in the x–y plane
   *
   *  Hence we are taking the triangles formed by (i, i+1, j+1)
   *  and (i, i-1, j-1) whose total area covers a unit cell,
   *  and ommit the other two triangles (bottom right and top left)
   *  When we sum this quantity across the whole lattice we cover every
   *  triangle. The final quantity will be scaled by 8 PI
   *  to match +-1 for a full skyrmion configuration
   *
   *  Recently, other ways to calculate a discrete skyrmion number have
   *  been proposed: http://arxiv.org/pdf/1601.08212.pdf
   *                 Phys. Rev. B 93, 024417
   *
   *  also based on using three spins using triangles. This could be
   *  useful for applying to a hexagonal lattice in the future.
   *
   */

  int i;
  int index, id;

  double sum = 0;

  double S[3], S_i[3], S_j[3];

  int nxy = nx * ny;

  for (i = 0; i < nxy; i++) {
    index = 3 * i;

    /* The starting index of the nearest neighbours for the
     * i-th spin */
    int id_nn = n_ngbs * i;

    S[0] = spin[index];
    S[1] = spin[index + 1];
    S[2] = spin[index + 2];

    S_i[0] = S_i[1] = S_i[2] = 0;
    S_j[0] = S_j[1] = S_j[2] = 0;

    // neighbour at -x
    // Remember that the index is -1 for sites without material
    if (ngbs[id_nn] > 0) {
      id = 3 * ngbs[id_nn];
      S_i[0] = spin[id];
      S_i[1] = spin[id + 1];
      S_i[2] = spin[id + 2];
    }

    // neighbour at -y
    if (ngbs[id_nn + 2] > 0) {
      id = 3 * ngbs[id_nn + 2];
      S_j[0] = spin[id];
      S_j[1] = spin[id + 1];
      S_j[2] = spin[id + 2];
    }

    // The  S_i \dot ( S_{i+1} \times S_{j+1} )
    charge[i] = volume(S, S_i, S_j);

    S_i[0] = S_i[1] = S_i[2] = 0;
    S_j[0] = S_j[1] = S_j[2] = 0;

    // neighbour at +x
    if (ngbs[id_nn + 1] > 0) {
      id = 3 * ngbs[id_nn + 1];
      S_i[0] = spin[id];
      S_i[1] = spin[id + 1];
      S_i[2] = spin[id + 2];
    }

    // neighbour at +y
    if (ngbs[id_nn + 3] > 0) {
      id = 3 * ngbs[id_nn + 3];
      S_j[0] = spin[id];
      S_j[1] = spin[id + 1];
      S_j[2] = spin[id + 2];
    }

    //  The S_i \dot ( S_{i-1} \times S_{j-1} )
    charge[i] += volume(S, S_i, S_j);

    /* Scale the chirality quantity */
    charge[i] /= (8 * WIDE_PI);

    /* We use the sum to output the total spin chirality
     * or skyrmion number */
    sum += charge[i];
  }

  return sum;
}

double compute_BergLuscher_angle(double *s1, double *s2, double *s3) {

  /* Compute the spherical angle given by the neighbouring spin 3-vectors s1,
   * s2 and s3 (these are arrays that can be also passed as pointers; see the
   * -dot- function). The spherical angle Omega, defined by the
   *  vectors in a unit sphere, is computed using the imaginary exponential
   *  defined by Berg and Luscher [Nucl Phys B 190, 412 (1981)]:
   *   exp (i sigma Omega) = 1 + s1 * s2 + s2 * s3 + s3 * s1 + i s1 * (s2 X s3)
   *                         ------------------------------------------------
   *                            2 (1 + s1 * s2) (1 + s2 * s3) (1 + s3 * s1)
   * The denominator is a normalisation factor which we call rho, and i
   * stands for an imaginary number in the numerator. The factor sigma is an
   * orientation given by sign(s1 * s2 X s3), however, we do not use it since
   * we take counter clock wise directions for the (s1, s2, s3) triangle of
   * spins, as pointed out by Yin et al. [PRB 93, 174403 (2016)].
   *
   * Therefore we use a complex logarithm:
   *      clog( r * exp(i theta) ) = r + i theta
   * to calculate the angle Omega, since the clog is well defined in the
   * [-PI, PI] range, giving the correct sign for the topological number (we
   * could also use the arcsin when decomposing the exp).
   *
   * Notice we normalise the angle by a 4 PI factor
   *
   */

  double rho;
  std::complex<double> exp;
  double crossp[3];

  cross(crossp, s2, s3);

  rho = sqrt(2 * (1 + dot(&s1[0], &s2[0], 3)) * (1 + dot(&s2[0], &s3[0], 3)) *
             (1 + dot(&s3[0], &s1[0], 3)));

  exp = (1 + dot(&s1[0], &s2[0], 3) + dot(&s2[0], &s3[0], 3) + dot(&s3[0], &s1[0], 3) +
         1i * dot(&s1[0], &crossp[0], 3)) /
        rho;

  return 2 * std::imag(std::log(exp)) / (4 * WIDE_PI);
}

double skyrmion_number_BergLuscher(double *spin, double *charge, int nx, int ny,
                                   int nz, int *ngbs, int n_ngbs) {

  /* Compute the topological charge (or skyrmion number) by adding triangles
   * of neighbouring spins for every lattice site, which cover triangle areas
   * in a unit sphere (i.e. we map the lattice area into a unit sphere
   * surface, using a triangulation). The exponential that defines every
   * spherical angle was firstly mentioned by Berg and Luscher [Nucl Phys B
   * 190, 412 (1981)] for a discrete square lattice but we generalise it here
   * for hexagonal crystals.
   *
   * NEIGHBOURS DEFINITION:
   *
   * *ngbs is the array with the neighbours information for every lattice
   * site. The first 6 elements are the NEAREST neighbours. For cuboid meshes,
   * the array is like:
   *
   *                                          __ indexes beyond nearest ngbs
   *                                          |
   *  NN:      j =0  j=1 ...                  |  j=0  j=1  ...
   *         [ 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, ... 1-x, 1+x, 1-y, ...  ]
   *  spin:   i=n_ngbs * 0                       i=n_ngbs * 1
   * ...
   *
   * where  0-y  is the index of the neighbour of the 0th spin, in the -y
   * direction, for example, so the neighbours of the i-th spin start at the
   * (n_ngbs * i) position of the array. This is similar for hexagonal meshes.
   *
   * For a cuboid mesh,  for the i-th spin, we generate the nearest ngbs
   * triangles using the triangles: [i, j=1, j=3] , [i, j=0, j=2]
   *                                     +x   +y         -x   -y
   * i.e. we cover the top right and bottom left areas ( we could also use
   * the top left and bottom right)
   *
   * For a hexagonal mesh,  for the i-th spin, we generate the nearest ngbs
   * triangles using the triangles: [i, j=0, j=2] , [i, j=0, j=2]
   *                                     E    NE         W    SW
   *                                   (East) ...
   * whose area covers a unit cell in a hexagonal crystal arrangement.
   *
   *
   * Since the NEAREST neighbours agree in indexes we can use the same function
   * for both meshes.
   *
   * ------------------------------------------------------------------------
   *
   * Having the triangles defined, we compute the spherical triangle area
   * spanned by the three spins using the compute_BergLuscher_angle function.
   *
   * The total topological charge Q is computed summing all triangles:
   *                 __
   *         Q  =   \     1   [ Omega(S_i, S_0, S_2) + Omega(S_i, S_1, S_3) ]
   *                /__  ---
   *                    4 PI
   *                 i
   *
   * ehich are normalised by 4 PI, so we get an integer number for the number
   * of times the unit sphere is covered by the spin directions in every
   * triangle, e.g. if we have two skyrmions we get approximately Q = 2.
   *
   */

  int n = nx * ny * nz;
  int i, spin_index;
  double total_sum = 0;

  // Sweep through every lattice site
  for (i = 0; i < n; i++) {

    // x-y-z components for the i-th spin start at:
    spin_index = 3 * i;

    // Reset the charge array
    charge[i] = 0;

    // Compute the spherical triangle area for the triangle formed
    // by the i-th spin and neighbours 0 and 2, i.e.
    // [i j=0 j=2]. First check that the NNs exist:
    if (ngbs[n_ngbs * i] >= 0 && ngbs[n_ngbs * i + 2] >= 0) {
      charge[i] +=
          compute_BergLuscher_angle(&spin[spin_index], &spin[3 * ngbs[n_ngbs * i]],
                                    &spin[3 * ngbs[n_ngbs * i + 2]]);
    }

    // Triangle: [i j=1 j=3]
    if (ngbs[n_ngbs * i + 1] >= 0 && ngbs[n_ngbs * i + 3] >= 0) {
      charge[i] += compute_BergLuscher_angle(&spin[spin_index],
                                             &spin[3 * ngbs[n_ngbs * i + 1]],
                                             &spin[3 * ngbs[n_ngbs * i + 3]]);
    }

    total_sum += charge[i];
  }

  return total_sum;
}

// compute the first derivative respect to x and for the whole mesh
// assume 2d pbc is used
void compute_px_py_c(double *spin, int nx, int ny, int nz, double *px,
                     double *py) {
  int nyz = ny * nz;
  int n1 = nx * nyz, n2 = 2 * n1;

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      int index = nyz * i + nz * j;

      // x-1
      int id1 = index - nyz;
      if (i == 0) {
        id1 += n1;
      }

      // x+1
      int id2 = index + nyz;
      if (i == nx - 1) {
        id2 -= n1;
      }

      px[index] = (spin[id2] - spin[id1]) / 2.0;
      px[index + n1] = (spin[id2 + n1] - spin[id1 + n1]) / 2.0;
      px[index + n2] = (spin[id2 + n2] - spin[id1 + n2]) / 2.0;

      // y-1
      id1 = index - nz;
      if (j == 0) {
        id1 += nyz;
      }

      // y+1
      id2 = index + nz;
      if (j == ny - 1) {
        id2 -= nyz;
      }

      py[index] = (spin[id2] - spin[id1]) / 2.0;
      py[index + n1] = (spin[id2 + n1] - spin[id1 + n1]) / 2.0;
      py[index + n2] = (spin[id2 + n2] - spin[id1 + n2]) / 2.0;
    }
  }
}

inline int get_index(int i, int j, int k, int nx, int nxy) {

  return k * nxy + j * nx + i;
}

// compute the guiding centre, Dynamics of magnetic vortices,     N.
// Papanicolaou,
// T.N. Tomaras 360, 425-462, (1991)
void compute_guiding_center(double *spin, int nx, int ny, int nz, int nx_start,
                            int nx_stop, int ny_start, int ny_stop,
                            double *res) {

  int nxy = ny * nx;
  int i, j;
  int id, index;

  double charge;
  double sum = 0, Rx = 0, Ry = 0;

  double S[3], S_i[3], S_j[3];
  int k = 0;

  for (i = nx_start; i < nx_stop; i++) {
    for (j = ny_start; j < ny_stop; j++) {
      index = 3 * get_index(i, j, k, nx, nxy);
      S[0] = spin[index];
      S[1] = spin[index + 1];
      S[2] = spin[index + 2];

      S_i[0] = S_i[1] = S_i[2] = 0;
      S_j[0] = S_j[1] = S_j[2] = 0;
      if (j > 0) {
        id = 3 * get_index(i, j - 1, k, nx, nxy);
        S_j[0] = spin[id];
        S_j[1] = spin[id + 1];
        S_j[2] = spin[id + 2];
      }

      if (i > 0) {
        id = 3 * get_index(i - 1, j, k, nx, nxy);
        S_i[0] = spin[id];
        S_i[1] = spin[id + 1];
        S_i[2] = spin[id + 2];
      }

      charge = volume(S, S_i, S_j);
      sum += charge;
      Rx += i * charge;
      Ry += j * charge;

      S_i[0] = S_i[1] = S_i[2] = 0;
      S_j[0] = S_j[1] = S_j[2] = 0;
      if (i < nx - 1) {
        id = 3 * get_index(i + 1, j, k, nx, nxy);
        S_i[0] = spin[id];
        S_i[1] = spin[id + 1];
        S_i[2] = spin[id + 2];
      }

      if (j < ny - 1) {
        id = 3 * get_index(i, j + 1, k, nx, nxy);
        S_j[0] = spin[id];
        S_j[1] = spin[id + 1];
        S_j[2] = spin[id + 2];
      }

      charge = volume(S, S_i, S_j);
      sum += charge;
      Rx += i * charge;
      Ry += j * charge;
    }
  }

  res[0] = Rx / sum;
  res[1] = Ry / sum;
}
