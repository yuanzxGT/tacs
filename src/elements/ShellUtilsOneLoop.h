/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_SHELL_UTILS_OneLoop_H
#define TACS_SHELL_UTILS_OneLoop_H

/*
  Common shell related computations.

  To evaluate the strain at the shell surface:

  1. Comute the transformation required at the current point
  2. Evaluate the strain/bmat/product of stress with
  the second derivative of the strain etc
*/

#include "ShellUtils.h"


/*
  The 'natural' shell transformations. These take the first local
  coordinate axis to be the line formed by the parametrization of the
  element surface. This is always defined for non-degenerate surfaces
*/

TACS_BEGIN_NAMESPACE(shellutils)




/*
  Compute the displacement and the derivative of the displacement
  along the parametric directions

  input:
  num_nodes: the number of nodes
  vars: the x,y,z displacements and rotations for each node
  N: the shape functions
  Na: the derivative of the shape functions along the first
  parametric direction
  Nb: the derivative of the shape functions along the second
  parametric direction

  output:
  U: the displacements and rotations at the parametric point
  Ud: the derivative of the displacements along the parametric directions
*/
void compute_shell_Ud_tri( const int num_nodes,
                       TacsScalar U[], TacsScalar Ud[],
                       const TacsScalar vars[],
                       const double N[], const double Na[],
                       const double Nb[] );


/*
  Compute the position and first derivatives along the shell surface.
  At the same time, compute the shape functions.

  The shell surface is based on the shape functions and nodal
  locations.

  input:
  gpt: the current Gauss point
  Xpts: the nodal locations for the

  output:
  X, Xd: position and first derivative of the surface
  N: the shape functions
  Na, Nb: the first derivatives of the shape functions
*/
void shell_jacobian_tri( const int order, TacsScalar X[],
                     TacsScalar Xd[],
                     double N[], double Na[], double Nb[],
                     const double gpt[],
                     const TacsScalar Xpts[] );


/*
  Compute the postion as well as the first and second derivatives
  along the surface of the shell.

  At the same time compute the shape functions, and their derivatives,
  at the current Gauss point gpt. Here Xpts is a vector of the nodal
  locations.

  input:
  gpt: the current Gauss point
  Xpts: the nodal locations for the

  output:
  X, Xd, Xdd: position, first and second derivatives of the surface
  N: the shape functions
  Na, Nb: the first derivatives of the shape functions
  Naa, Nab, Nbb: the second derivatives of the shape functions
*/
void shell_hessian_tri( const int order,
                    TacsScalar X[], TacsScalar Xd[],
                    TacsScalar Xdd[],
                    double N[], double Na[], double Nb[],
                    double Naa[], double Nab[], double Nbb[],
                    const double gpt[],
                    const TacsScalar Xpts[] );

/*
  Compute only the terms/derivatives required for computing the
  in-plane tensorial strain components.

  input:
  order:     order of the interpolation
  pt:        the parametric point to evaluate the shape functions
  Xpts:      the nodal locations
  vars:      the state varaible values

  output:
  Urot: the rotations at the mid-surface
  Ud: the derivatives of the displacements at the mid-surface
  Xd: the derivatives of the nodal locations along the parametric directions
  N, Na, Nb:  the tensor product shape functions and their derivatives
*/
void compute_tensorial_components_tri( const int order,
                                   TacsScalar Urot[],
                                   TacsScalar Ud[], TacsScalar Xd[],
                                   double N[], double Na[], double Nb[],
                                   const double pt[],
                                   const TacsScalar Xpts[],
                                   const TacsScalar vars[] );


















/*
  Compute the second derivative of the strain w.r.t. the
  displacements/rotations - many of these components are the same.

  input:
  N11, N22, N12: the shape functions associated with the interpolation
  of the strain components

  knots, pknots: the tying interpolation points

  ouput:
  n11, n22, n12: the second derivative of the strain w.r.t. the nodal
  displacements
*/
template <int order, int tying_order>
void nonlinear_tying_nmat_tri( TacsScalar n11[], TacsScalar n22[],
                           TacsScalar n12[],
                           const double N11[], const double N22[],
                           const double N12[],
                           const double knots[][2], const double pknots[][2] ){
  double N[3], Na[3], Nb[3];

  // Evaluate g11 g22 g12 at GaussPts
  for ( int m = 0; m < tying_order; m++ ){
      double pt[2];
      pt[0] = pknots[m][0];
      pt[1] = pknots[m][1];

      N[0] = 1.0 - (pt[0] + pt[1]);
	  N[1] = pt[0];
	  N[2] = pt[1];

      // First derivatives
      Na[0] = - 1.0;
      Na[1] = 1.0;
      Na[2] = 0.0;

      Nb[0] = - 1.0;
      Nb[1] = 0.0;
      Nb[2] = 1.0;

      TacsScalar * _n22 = n22;
	  TacsScalar * _n12 = n12;
	  TacsScalar * _n11 = n11;
      for ( int i = 0; i < 3; i++ ){
        for ( int j = 0; j <= i; j++ ){
          _n22[0] += Nb[i]*Nb[j]*N22[0];
          _n22++;
		  
		  _n12[0] += 0.5*(Na[i]*Nb[j] + Nb[i]*Na[j])*N12[0];
          _n12++;
		  
		  _n11[0] += Na[i]*Na[j]*N11[0];
          _n11++;
        }
      }
      N22++;
	  N12++;
	  N11++;
    
  }

}

/*
  Add the derivative of the product of stress and the derivative of
  the strain to part of the matrix.  This is used for the nonlinear
  contributions from the in-plane components of the tensorial strain.
*/
template <int order, int tying_order>
void add_nonlinear_tying_stress_nmat_tri( TacsScalar matrix[],
                                      TacsScalar scale,
                                      const TacsScalar stress[],
                                      const TacsScalar tx[],
                                      const double N11[],
                                      const double N22[],
                                      const double N12[],
                                      const double knots[][2],
                                      const double pknots[][2],
                                      const TacsScalar Xpts[] ){
  const int nvars = 6*order;
  const int size = (order*(order+1))/2;

  TacsScalar n11[size], n22[size], n12[size];
  memset(n11, 0, size*sizeof(TacsScalar));
  memset(n22, 0, size*sizeof(TacsScalar));
  memset(n12, 0, size*sizeof(TacsScalar));

  nonlinear_tying_nmat_tri<order, tying_order>(n11, n22, n12,
                                           N11, N22, N12, knots, pknots);

  TacsScalar * _n11 = n11;
  TacsScalar * _n12 = n12;
  TacsScalar * _n22 = n22;

  // Compute the second derivative contributions
  for ( int i = 0; i < order; i++ ){
    const int row = 6*i;
    for ( int j = 0; j <= i; j++ ){
      const int col = 6*j;

      // Evaluate the second derivative of the tensorial strain.
      // Note that this is NOT the engineering strain!
      TacsScalar g[6], s[6];
      g[0] = _n11[0];
      g[1] = _n22[0];
      g[5] = _n12[0];
      g[2] = g[3] = g[4] = 0.0;

      // Transform to the local axis - note that this employs the
      // transformation for stress - since we're using tensorial strain -
      // not engineering strain
      Tensor::transform3DStress(s, g, tx);

      // Evaluate the contribution - note the factor of 2.0 due to the
      // use of tensorial strain
      TacsScalar val =
        scale*(stress[0]*s[0] + stress[1]*s[1] +
               2.0*(stress[2]*s[5] + stress[6]*s[3] + stress[7]*s[4]));

      // Add values to the diagonal of each component
      for ( int k = 0; k < 3; k++ ){
        matrix[(row+k)*nvars + col+k] += val;
      }

      // Increment the pointers
      _n11++; _n22++; _n12++;
    }
  }
}





/*
  Evaluate the tensorial displacement-based strain at the tying points.

  The tensorial strain is given as follows,

  E_t = U_{x,xi}*X_{,xi}^{T} + X_{,xi}*U_{x,xi}^{T}

  X_{,xi} = Xd^{T}
  U_{x,xi} = [ Ud, r ]

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  vars: the displacements and rotations at the nodal locations
  Xpts: the nodal locations for the element

  output:
  g11: the in-plane normal tensorial strain in the 1 direction
  g22: the in-plane normal tensorial strain in the 2 direction
  g12: the in-plane tensorial shear strain

  g23: the out-of-plane tensorial shear strain in the 2-direction
  g13: the out-of-plane tensorial shear strain in the 1-direction
*/
template <int order, int tying_order>
void compute_tying_strain_tri( const int is_linear,
                           TacsScalar g11[], TacsScalar g22[],
                           TacsScalar g12[],
                           TacsScalar g23[], TacsScalar g13[],
                           const double knots[][2], const double pknots[][2],
                           const TacsScalar vars[],
                           const TacsScalar Xpts[] ){
  double N[3], Na[3], Nb[3];

  // Evaluate g13 and g23
  // MITC 3 only mixed interpolation for transverse shear g23 g13 at tying point
  for ( int m = 0; m < tying_order; m++ ){
      double pt[2];
      pt[0] = knots[m][0];
      pt[1] = knots[m][1];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components_tri(order, Urot, Ud, Xd,
                                   N, Na, Nb, pt, Xpts, vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation
      TacsScalar r[3];
      Tensor::crossProduct3D(r, Urot, normal);
      
      if (is_linear){
		  

        g23[0] = 0.5*(Xd[3]*r[0] + Xd[4]*r[1] + Xd[5]*r[2] +
                      normal[0]*Ud[1] + normal[1]*Ud[3] + normal[2]*Ud[5]);
		g13[0] = 0.5*(Xd[0]*r[0] + Xd[1]*r[1] + Xd[2]*r[2] +
                      normal[0]*Ud[0] + normal[1]*Ud[2] + normal[2]*Ud[4]);
         g23++; g13++;

      }
      else {
        g23[0] = 0.5*(Xd[3]*r[0] + Xd[4]*r[1] + Xd[5]*r[2] +
                      normal[0]*Ud[1] + normal[1]*Ud[3] + normal[2]*Ud[5]);
		g13[0] = 0.5*(Xd[0]*r[0] + Xd[1]*r[1] + Xd[2]*r[2] +
                      normal[0]*Ud[0] + normal[1]*Ud[2] + normal[2]*Ud[4]);
        g23++; g13++;

      }
    
  }

  // Evaluate in plane strain g11, g22, g12, at GaussPts
  for ( int m = 0; m < tying_order; m++ ){
      double pt[2];
      pt[0] = pknots[m][0];
      pt[1] = pknots[m][1];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components_tri(order, Urot, Ud, Xd,
                                   N, Na, Nb, pt, Xpts, vars);

      if (is_linear){
        g12[0] = 0.5*(Xd[0]*Ud[1] + Xd[1]*Ud[3] + Xd[2]*Ud[5] +
                      Xd[3]*Ud[0] + Xd[4]*Ud[2] + Xd[5]*Ud[4]);
		g22[0] = Xd[3]*Ud[1] + Xd[4]*Ud[3] + Xd[5]*Ud[5];			  
		g11[0] = Xd[0]*Ud[0] + Xd[1]*Ud[2] + Xd[2]*Ud[4];
        g11++;
        g12++; g22++;

      }
      else {
        g12[0] = 0.5*(Xd[0]*Ud[1] + Xd[1]*Ud[3] + Xd[2]*Ud[5] +
                      Xd[3]*Ud[0] + Xd[4]*Ud[2] + Xd[5]*Ud[4] +
                      Ud[0]*Ud[1] + Ud[2]*Ud[3] + Ud[4]*Ud[5]);
		g22[0] = (Xd[3]*Ud[1] + Xd[4]*Ud[3] + Xd[5]*Ud[5] +
                  0.5*(Ud[1]*Ud[1] + Ud[3]*Ud[3] + Ud[5]*Ud[5]));
		g11[0] = (Xd[0]*Ud[0] + Xd[1]*Ud[2] + Xd[2]*Ud[4] +
                  0.5*(Ud[0]*Ud[0] + Ud[2]*Ud[2] + Ud[4]*Ud[4]));     
        g12++; g22++; g11++;


      }
    
  }
}





/*
  Compute the linear tying strain bmat matrices. This corresponds to
  the derivative of the displacement-based strain at each of the tying
  points.

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  Xpts: the nodal locations for the element

  output:
  b11: the derivative of the in-plane normal tensorial strain in the 1
  direction
  b22: the derivative of the in-plane normal tensorial strain in the 2
  direction
  b12: the derivative of the in-plane tensorial shear strain

  b23: the derivative of the out-of-plane tensorial shear strain in
  the 2-direction
  b13: the derivative of the out-of-plane tensorial shear strain in
  the 1-direction
*/
template <int order, int tying_order>
void compute_tying_bmat_tri( const int is_linear,
                         TacsScalar g11[], TacsScalar g22[],
                         TacsScalar g12[],
                         TacsScalar g23[], TacsScalar g13[],
                         TacsScalar b11[], TacsScalar b22[],
                         TacsScalar b12[],
                         TacsScalar b23[], TacsScalar b13[],
                         const double knots[][2], const double pknots[][2],
                         const TacsScalar vars[],
                         const TacsScalar Xpts[] ){
  double N[3], Na[3], Nb[3];

  // Evaluate g13 and g23
  // MITC 3 only mixed interpolation for transverse shear g23 g13 at tying point
  for ( int m = 0; m < 3; m++ ){
      double pt[2];
      pt[0] = knots[m][0];
      pt[1] = knots[m][1];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components_tri(order, Urot, Ud, Xd,
                                   N, Na, Nb, pt, Xpts, vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation
      TacsScalar r[3];
      Tensor::crossProduct3D(r, Urot, normal);
      
      if (is_linear){
		  

        g23[0] = 0.5*(Xd[3]*r[0] + Xd[4]*r[1] + Xd[5]*r[2] +
                      normal[0]*Ud[1] + normal[1]*Ud[3] + normal[2]*Ud[5]);
		g13[0] = 0.5*(Xd[0]*r[0] + Xd[1]*r[1] + Xd[2]*r[2] +
                      normal[0]*Ud[0] + normal[1]*Ud[2] + normal[2]*Ud[4]);
         g23++; g13++;

        // Add the derivative to the b-matrix
        for ( int i = 0; i < order; i++ ){


          b23[0] = 0.5*normal[0]*Nb[i];
          b23[1] = 0.5*normal[1]*Nb[i];
          b23[2] = 0.5*normal[2]*Nb[i];
          b23[3] = 0.5*N[i]*(Xd[5]*normal[1] - Xd[4]*normal[2]);
          b23[4] = 0.5*N[i]*(Xd[3]*normal[2] - Xd[5]*normal[0]);
          b23[5] = 0.5*N[i]*(Xd[4]*normal[0] - Xd[3]*normal[1]);
          b23 += 6;
		  
		  b13[0] = 0.5*normal[0]*Na[i];
          b13[1] = 0.5*normal[1]*Na[i];
          b13[2] = 0.5*normal[2]*Na[i];
          b13[3] = 0.5*N[i]*(Xd[2]*normal[1] - Xd[1]*normal[2]);
          b13[4] = 0.5*N[i]*(Xd[0]*normal[2] - Xd[2]*normal[0]);
          b13[5] = 0.5*N[i]*(Xd[1]*normal[0] - Xd[0]*normal[1]);
          b13 += 6;
        }
      }
      else {
        g23[0] = 0.5*(Xd[3]*r[0] + Xd[4]*r[1] + Xd[5]*r[2] +
                      normal[0]*Ud[1] + normal[1]*Ud[3] + normal[2]*Ud[5]);
		g13[0] = 0.5*(Xd[0]*r[0] + Xd[1]*r[1] + Xd[2]*r[2] +
                      normal[0]*Ud[0] + normal[1]*Ud[2] + normal[2]*Ud[4]);
        g23++; g13++;

        // Add the derivative to the b-matrix
        for ( int i = 0; i < order; i++ ){

          b23[0] = 0.5*normal[0]*Nb[i];
          b23[1] = 0.5*normal[1]*Nb[i];
          b23[2] = 0.5*normal[2]*Nb[i];
          b23[3] = 0.5*N[i]*(Xd[5]*normal[1] - Xd[4]*normal[2]);
          b23[4] = 0.5*N[i]*(Xd[3]*normal[2] - Xd[5]*normal[0]);
          b23[5] = 0.5*N[i]*(Xd[4]*normal[0] - Xd[3]*normal[1]);
          b23 += 6;
		  
		  b13[0] = 0.5*normal[0]*Na[i];
          b13[1] = 0.5*normal[1]*Na[i];
          b13[2] = 0.5*normal[2]*Na[i];
          b13[3] = 0.5*N[i]*(Xd[2]*normal[1] - Xd[1]*normal[2]);
          b13[4] = 0.5*N[i]*(Xd[0]*normal[2] - Xd[2]*normal[0]);
          b13[5] = 0.5*N[i]*(Xd[1]*normal[0] - Xd[0]*normal[1]);
          b13 += 6;
        }
      }
    
  }

  // Evaluate in plane strain g11, g22, g12, at GaussPts
  for ( int m = 0; m < 3; m++ ){
      double pt[2];
      pt[0] = pknots[m][0];
      pt[1] = pknots[m][1];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components_tri(order, Urot, Ud, Xd,
                                   N, Na, Nb, pt, Xpts, vars);

      if (is_linear){
        g12[0] = 0.5*(Xd[0]*Ud[1] + Xd[1]*Ud[3] + Xd[2]*Ud[5] +
                      Xd[3]*Ud[0] + Xd[4]*Ud[2] + Xd[5]*Ud[4]);
		g22[0] = Xd[3]*Ud[1] + Xd[4]*Ud[3] + Xd[5]*Ud[5];			  
		g11[0] = Xd[0]*Ud[0] + Xd[1]*Ud[2] + Xd[2]*Ud[4];
        g11++;
        g12++; g22++;

        for ( int i = 0; i < 3; i++ ){
          b12[0] = 0.5*(Xd[3]*Na[i] + Xd[0]*Nb[i]);
          b12[1] = 0.5*(Xd[4]*Na[i] + Xd[1]*Nb[i]);
          b12[2] = 0.5*(Xd[5]*Na[i] + Xd[2]*Nb[i]);
          b12 += 3;
		  
		  b22[0] = Xd[3]*Nb[i];
          b22[1] = Xd[4]*Nb[i];
          b22[2] = Xd[5]*Nb[i];
          b22 += 3;
		  
		  b11[0] = Xd[0]*Na[i];
          b11[1] = Xd[1]*Na[i];
          b11[2] = Xd[2]*Na[i];
          b11 += 3;
        }
      }
      else {
        g12[0] = 0.5*(Xd[0]*Ud[1] + Xd[1]*Ud[3] + Xd[2]*Ud[5] +
                      Xd[3]*Ud[0] + Xd[4]*Ud[2] + Xd[5]*Ud[4] +
                      Ud[0]*Ud[1] + Ud[2]*Ud[3] + Ud[4]*Ud[5]);
		g22[0] = (Xd[3]*Ud[1] + Xd[4]*Ud[3] + Xd[5]*Ud[5] +
                  0.5*(Ud[1]*Ud[1] + Ud[3]*Ud[3] + Ud[5]*Ud[5]));
		g11[0] = (Xd[0]*Ud[0] + Xd[1]*Ud[2] + Xd[2]*Ud[4] +
                  0.5*(Ud[0]*Ud[0] + Ud[2]*Ud[2] + Ud[4]*Ud[4]));     
        g12++; g22++; g11++;

        for ( int i = 0; i < 3; i++ ){
          b12[0] = 0.5*((Xd[3] + Ud[1])*Na[i] + (Xd[0] + Ud[0])*Nb[i]);
          b12[1] = 0.5*((Xd[4] + Ud[3])*Na[i] + (Xd[1] + Ud[2])*Nb[i]);
          b12[2] = 0.5*((Xd[5] + Ud[5])*Na[i] + (Xd[2] + Ud[4])*Nb[i]);
          b12 += 3;
		  
		  b22[0] = (Xd[3] + Ud[1])*Nb[i];
          b22[1] = (Xd[4] + Ud[3])*Nb[i];
          b22[2] = (Xd[5] + Ud[5])*Nb[i];
          b22 += 3;
		  
		  b11[0] = (Xd[0] + Ud[0])*Na[i];
          b11[1] = (Xd[1] + Ud[2])*Na[i];
          b11[2] = (Xd[2] + Ud[4])*Na[i];
          b11 += 3;
        }
      }
    
  }


}




/*
  Add the tying strain contribution to the strains. This adds the
  result of the interpolation of the displacement-based strain from
  the tying points to the strain. Only the in-plane and shear strain
  components are added.

  input:
  tx: the transformation to the local coordinates
  g11, g22, g12, g23, g13: the components of the tying strain

  N11, N22, N12: the interpolations

  output:
  strain: the strain values are set into the in-plane components (0, 1, 2)
  and the shear strain components (6, 7)
*/
template <int tying_order>
void add_tying_strain_tri( TacsScalar strain[], const double pt[], const TacsScalar tx[],
                       const TacsScalar g11[], const TacsScalar g22[],
                       const TacsScalar g12[],
                       const TacsScalar g23[], const TacsScalar g13[],
                       const double N11[], const double N22[],
                       const double N12[] ){
  TacsScalar g[6] = {0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0};

  // Use the interpolations
  
  
  
  for ( int k = 0; k < 3; k++ ){
    g[5] += g12[k]*N12[k];
    g[0] += g11[k]*N11[k];
    g[1] += g22[k]*N22[k];
  }  
  
  // transverse shear
  
  double cc;
  cc = g23[1] - g13[0] - g23[2] + g13[2];
  g[4] = g13[0] + cc * pt[1];
  g[3] = g23[1] - cc * pt[0];
  




  // Transform the tensorial strain g, to the local axis
  TacsScalar s[6];
  Tensor::transform3DStress(s, g, tx);
  strain[0] = s[0];
  strain[1] = s[1];
  strain[2] = 2.0*s[5];
  strain[6] = 2.0*s[3];
  strain[7] = 2.0*s[4];
}


/*
  Add the B matrices (the derivative of the strain w.r.t. the nodal
  displacements and rotations at each tying point) together weighted
  by the tying interpolations: N11, N22, N22. Transform the tensorial
  strain into the local strain.

  input:
  tx: the transformation to the local coordinates
  b11, b22, b12, b23, b13: the tying strain matrices
  N11, N22, N12: the tying strain shape functions

  output:
  B: the derivative of the strain w.r.t. the nodal displacements - this
  only writes in the in-plane and out-of-plane shear components (not
  the bending components.)
*/
template <int tying_order>
void add_tying_bmat_tri( TacsScalar B[], const double pt[], const int num_nodes,
                     const TacsScalar tx[],
                     const TacsScalar b11[], const TacsScalar b22[],
                     const TacsScalar b12[],
                     const TacsScalar b23[], const TacsScalar b13[],
                     const double N11[], const double N22[],
                     const double N12[] ){
  // Variables stored by b11[nvars*(n + m*(order-1)) + row] Points stored
  const int n6pts = 6*num_nodes;
  const int n3pts = 3*num_nodes;


     double cc;


  for ( int i = 0; i < num_nodes; i++ ){
    for ( int ii = 0; ii < 6; ii++ ){
      TacsScalar g[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

      int row = 6*i + ii;
      if ( ii < 3 ){
        int row3 = 3*i + ii;
        // at k-th tying point
        for ( int k = 0; k < 3; k++ ){
		  g[5] += b12[n3pts*k + row3]*N12[k];
          g[0] += b11[n3pts*k + row3]*N11[k];
          g[1] += b22[n3pts*k + row3]*N22[k];
		}  

		 cc=b23[n6pts*1 + row] - b13[n6pts*0 + row] - b23[n6pts*2 + row] + b13[n6pts*2 + row];
		 g[4] += b13[n6pts*0 + row] + cc * pt[1];
		 g[3] += b23[n6pts*1 + row] - cc * pt[0];
		  
        
      }
      else {
		 cc=b23[n6pts*1 + row] - b13[n6pts*0 + row] - b23[n6pts*2 + row] + b13[n6pts*2 + row];
		 g[4] += b13[n6pts*0 + row] + cc * pt[1];
		 g[3] += b23[n6pts*1 + row] - cc * pt[0];
      }

      // Transform the tensorial strain g, to the local axis
      TacsScalar s[6];
      Tensor::transform3DStress(s, g, tx);
      B[0] = s[0];
      B[1] = s[1];
      B[2] = 2.0*s[5];
      B[6] = 2.0*s[3];
      B[7] = 2.0*s[4];
      B += 8;
    }
  }
}


/*
  Evaluate the tying-interpolation at the appropriate points
*/
template <int tying_order>
void tying_interpolation_tri( const double pt[],
                          double N11[], double N22[], double N12[],
                          const double knots[][2], const double pknots[][2] ){
  double na[tying_order], nb[tying_order];
  
  //  
  
  double r1, s1, r2, s3, r,s;
  
  r = pt[0];
  s = pt[1];
  r1 = pknots[0][0];
  s1 = pknots[0][1];
  r2 = pknots[1][0];
  s3 = pknots[2][1];
  


  // in-plane: no interpolation, use strains ar GaussPts directly
  // g11 g22 g12
  /*double errTol=0.00001;

  N12[0] = N12[1] = N12[2]=0;
  for( int i=0; i<3; i++ ){
     if ( abs( pt[0] - pknots[i][0] ) < errTol && abs( pt[1] - pknots[i][1] ) < errTol ){
      N12[i] =  1.0;
	  N22[i] =  1.0;
	  N11[i] =  1.0;
	}
  }
    
*/

   N12[0] = 1 - ( r - r1 ) / ( r2 - r1) - ( s - s1 ) / ( s3 - s1);
   N12[1] = ( r - r1 ) / ( r2 - r1);
   N12[2] = ( s - s1 ) / ( s3 - s1);
   N22[0]=N11[0]=N12[0];
   N22[1]=N11[1]=N12[1];
   N22[2]=N11[2]=N12[2];
}


TACS_END_NAMESPACE

#endif
