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

#include "FElibrary.h"
#include "ShellUtilsOneLoop.h"

/*
  Common structural shell utilities for small strains (Triangular Element with one Loop)
*/

TACS_BEGIN_NAMESPACE(shellutils)


/*
  This is a helper function for the large rotation strain utilities

  Given the order and the tensor-product contributions to the strain
  moments, compute the displacements, rotations and rate of change of
  displacements.
*/
void compute_tensorial_components_tri( const int order,
                                   TacsScalar Urot[], TacsScalar Ud[], TacsScalar Xd[],
                                   double N[], double Na[], double Nb[], const double pt[],
                                   const TacsScalar Xpts[], const TacsScalar vars[] ){
  if (order == 3){
      
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


    // Compute the rotations
    Urot[0] = vars[3]*N[0] + vars[9]*N[1]  + vars[15]*N[2];
    Urot[1] = vars[4]*N[0] + vars[10]*N[1] + vars[16]*N[2];
    Urot[2] = vars[5]*N[0] + vars[11]*N[1] + vars[17]*N[2];

    // Compute the derivatives of the displacement
    Ud[0]  = vars[0]*Na[0] + vars[6]*Na[1]  + vars[12]*Na[2];
    Ud[2]  = vars[1]*Na[0] + vars[7]*Na[1]  + vars[13]*Na[2];
    Ud[4]  = vars[2]*Na[0] + vars[8]*Na[1]  + vars[14]*Na[2];

    // Compute the derivative of X along a
    Xd[0] = Xpts[0]*Na[0] + Xpts[3]*Na[1] + Xpts[6]*Na[2];
    Xd[1] = Xpts[1]*Na[0] + Xpts[4]*Na[1] + Xpts[7]*Na[2];
    Xd[2] = Xpts[2]*Na[0] + Xpts[5]*Na[1] + Xpts[8]*Na[2];

    // Compute the derivative of the displacement
    Ud[1]  = vars[0]*Nb[0] + vars[6]*Nb[1]  + vars[12]*Nb[2];
    Ud[3]  = vars[1]*Nb[0] + vars[7]*Nb[1]  + vars[13]*Nb[2];
    Ud[5]  = vars[2]*Nb[0] + vars[8]*Nb[1]  + vars[14]*Nb[2];

    // Compute the derivative of X along b
    Xd[3] = Xpts[0]*Nb[0] + Xpts[3]*Nb[1] + Xpts[6]*Nb[2];
    Xd[4] = Xpts[1]*Nb[0] + Xpts[4]*Nb[1] + Xpts[7]*Nb[2];
    Xd[5] = Xpts[2]*Nb[0] + Xpts[5]*Nb[1] + Xpts[8]*Nb[2];
  }



}










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
                       const double Nb[] ){
  if (num_nodes == 3){
    U[0] = vars[0]*N[0] + vars[6]*N[1]  + vars[12]*N[2];
    U[1] = vars[1]*N[0] + vars[7]*N[1]  + vars[13]*N[2];
    U[2] = vars[2]*N[0] + vars[8]*N[1]  + vars[14]*N[2];
    U[3] = vars[3]*N[0] + vars[9]*N[1]  + vars[15]*N[2];
    U[4] = vars[4]*N[0] + vars[10]*N[1] + vars[16]*N[2];
    U[5] = vars[5]*N[0] + vars[11]*N[1] + vars[17]*N[2];

    Ud[0]  = vars[0]*Na[0] + vars[6]*Na[1]  + vars[12]*Na[2];
    Ud[2]  = vars[1]*Na[0] + vars[7]*Na[1]  + vars[13]*Na[2];
    Ud[4]  = vars[2]*Na[0] + vars[8]*Na[1]  + vars[14]*Na[2];
    Ud[6]  = vars[3]*Na[0] + vars[9]*Na[1]  + vars[15]*Na[2];
    Ud[8]  = vars[4]*Na[0] + vars[10]*Na[1] + vars[16]*Na[2];
    Ud[10] = vars[5]*Na[0] + vars[11]*Na[1] + vars[17]*Na[2];

    Ud[1]  = vars[0]*Nb[0] + vars[6]*Nb[1]  + vars[12]*Nb[2];
    Ud[3]  = vars[1]*Nb[0] + vars[7]*Nb[1]  + vars[13]*Nb[2];
    Ud[5]  = vars[2]*Nb[0] + vars[8]*Nb[1]  + vars[14]*Nb[2];
    Ud[7]  = vars[3]*Nb[0] + vars[9]*Nb[1]  + vars[15]*Nb[2];
    Ud[9]  = vars[4]*Nb[0] + vars[10]*Nb[1] + vars[16]*Nb[2];
    Ud[11] = vars[5]*Nb[0] + vars[11]*Nb[1] + vars[17]*Nb[2];
  }
 
}


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
void shell_jacobian_tri( const int order,
                     TacsScalar X[], TacsScalar Xd[],
                     double N[], double Na[], double Nb[],
                     const double gpt[],
                     const TacsScalar Xpts[] ){
  double na[8], nb[8];
  double dna[8], dnb[8];

  X[0] = X[1] = X[2] = TacsScalar(0.0);
  Xd[0] = Xd[1] = Xd[2] = TacsScalar(0.0);
  Xd[3] = Xd[4] = Xd[5] = TacsScalar(0.0);

  if (order <= 8){
    if (order == 3){
	  N[0] = 1.0 - (gpt[0] + gpt[1]);
	  N[1] = gpt[0];
	  N[2] = gpt[1];

      // First derivatives
      Na[0] = - 1.0;
      Na[1] = 1.0;
      Na[2] = 0.0;

      Nb[0] = - 1.0;
      Nb[1] = 0.0;
      Nb[2] = 1.0;


    }
    


    const int order2 = 3;
    for ( int i = 0; i < order2; i++ ){
      // Evaluate the derivatives
      X[0] += Xpts[0]*N[0];
      X[1] += Xpts[1]*N[0];
      X[2] += Xpts[2]*N[0];

      // First derivatives of A and B
      Xd[0] += Xpts[0]*Na[0];
      Xd[1] += Xpts[1]*Na[0];
      Xd[2] += Xpts[2]*Na[0];

      Xd[3] += Xpts[0]*Nb[0];
      Xd[4] += Xpts[1]*Nb[0];
      Xd[5] += Xpts[2]*Nb[0];

      N++; Na++; Nb++;
      Xpts += 3;
    }
  }
}




/*
  Compute the postion as well ass the first and second derivatives
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
                    const TacsScalar Xpts[] ){
  double na[8], nb[8];
  double dna[8], dnb[8];
  double ddna[8], ddnb[8];

  X[0] = X[1] = X[2] = TacsScalar(0.0);

  Xd[0] = Xd[1] = Xd[2] = TacsScalar(0.0);
  Xd[3] = Xd[4] = Xd[5] = TacsScalar(0.0);

  Xdd[0] = Xdd[1] = Xdd[2] = TacsScalar(0.0);
  Xdd[3] = Xdd[4] = Xdd[5] = TacsScalar(0.0);
  Xdd[6] = Xdd[7] = Xdd[8] = TacsScalar(0.0);

  if (order <= 8){
    if (order == 3){
	  N[0] = 1.0 - (gpt[0] + gpt[1]);
	  N[1] = gpt[0];
	  N[2] = gpt[1];

      // First derivatives
      Na[0] = - 1.0;
      Na[1] = 1.0;
      Na[2] = 0.0;

      Nb[0] = - 1.0;
      Nb[1] = 0.0;
      Nb[2] = 1.0;

      // Second derivatives
      Naa[0] = Naa[1] = Naa[2] = 0.0;
	  Nab[0] = Nab[1] = Nab[2] = 0.0;
      Nbb[0] = Nbb[1] = Nbb[2] = 0.0;
    }

    

    const int order2 = 3;
    for ( int i = 0; i < order2; i++ ){
      X[0] += Xpts[0]*N[0];
      X[1] += Xpts[1]*N[0];
      X[2] += Xpts[2]*N[0];

      // First derivatives
      Xd[0] += Xpts[0]*Na[0];
      Xd[1] += Xpts[1]*Na[0];
      Xd[2] += Xpts[2]*Na[0];

      Xd[3] += Xpts[0]*Nb[0];
      Xd[4] += Xpts[1]*Nb[0];
      Xd[5] += Xpts[2]*Nb[0];

      // Second derivatives
      Xdd[0] += Xpts[0]*Naa[0];
      Xdd[1] += Xpts[1]*Naa[0];
      Xdd[2] += Xpts[2]*Naa[0];

      Xdd[3] += Xpts[0]*Nab[0];
      Xdd[4] += Xpts[1]*Nab[0];
      Xdd[5] += Xpts[2]*Nab[0];

      Xdd[6] += Xpts[0]*Nbb[0];
      Xdd[7] += Xpts[1]*Nbb[0];
      Xdd[8] += Xpts[2]*Nbb[0];

      N++; Na++; Nb++;
      Naa++; Nab++; Nbb++;
      Xpts += 3;
    }
  }
}

TACS_END_NAMESPACE
