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

#include "TACSCompliance.h"
#include "TACSAssembler.h"

/*
  Initialize the Compliance class properties
*/
TACSCompliance::TACSCompliance( TACSAssembler *_assembler ):
  TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
               TACSFunction::SINGLE_STAGE, 0){
  compliance = 0.0;
}

TACSCompliance::~TACSCompliance(){}

const char *TACSCompliance::funcName = "TACSCompliance";

const char *TACSCompliance::getObjectName(){
  return funcName;
}

/*
  Retrieve the function value
*/
TacsScalar TACSCompliance::getFunctionValue(){
  return compliance;
}

/*
  Set the compliance to zero on all MPI processes
*/
void TACSCompliance::initEvaluation( EvaluationType ftype ){
  compliance = 0.0;
}

/*
  Sum the compliance across all MPI processors
*/
void TACSCompliance::finalEvaluation( EvaluationType ftype ){
  // Distribute the values of the KS function computed on this domain
  TacsScalar temp = compliance;
  MPI_Allreduce(&temp, &compliance, 1, TACS_MPI_TYPE,
                MPI_SUM, assembler->getMPIComm());
}

/*
  Evaluate the compliance contributed by this element
*/
void TACSCompliance::elementWiseEval( EvaluationType ftype,
                                      int elemIndex,
                                      TACSElement *element,
                                      double time,
                                      TacsScalar scale,
                                      const TacsScalar Xpts[],
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[],
                                      const TacsScalar ddvars[] ){
  // Retrieve the basis object for this element (if it exists)
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      // Evaluate the strain energy density
      TacsScalar U0 = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_STRAIN_ENERGY_DENSITY,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars, &U0);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        compliance += detJ*weight*U0;
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function to the state variables.
*/
void TACSCompliance::getElementSVSens( int elemIndex, TACSElement *element,
                                       double time,
                                       TacsScalar alpha,
                                       TacsScalar beta,
                                       TacsScalar gamma,
                                       const TacsScalar Xpts[],
                                       const TacsScalar vars[],
                                       const TacsScalar dvars[],
                                       const TacsScalar ddvars[],
                                       TacsScalar dfdu[] ){
  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->getNumVariables();
  memset(dfdu, 0, numVars*sizeof(TacsScalar));

  // Retrieve the basis object for this element (if it exists)
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      // Evaluate the strain energy density
      TacsScalar U0 = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_STRAIN_ENERGY_DENSITY,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars, &U0);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Evaluate the derivative of the strain energy w.r.t. state variables
        TacsScalar dfdq = detJ*weight;
        element->addPointQuantitySVSens(elemIndex, TACS_STRAIN_ENERGY_DENSITY, time,
                                        alpha, beta, gamma,
                                        i, pt, Xpts, vars, dvars, ddvars,
                                        &dfdq, dfdu);
      }
    }
  }
}

/*
  Retrieve the element contribution to the derivative of the function
  w.r.t. the element nodes
*/
void TACSCompliance::getElementXptSens( int elemIndex,
                                        TACSElement *element,
                                        double time,
                                        TacsScalar scale,
                                        const TacsScalar Xpts[],
                                        const TacsScalar vars[],
                                        const TacsScalar dvars[],
                                        const TacsScalar ddvars[],
                                        TacsScalar dfdXpts[] ){
  // Zero the derivative of the function w.r.t. the element node
  // locations
  int numNodes = element->getNumNodes();
  memset(dfdXpts, 0, 3*numNodes*sizeof(TacsScalar));
}

/*
  Evaluate the derivative of the compliance w.r.t. the material
  design variables
*/
void TACSCompliance::addElementDVSens( int elemIndex,
                                       TACSElement *element,
                                       double time,
                                       TacsScalar scale,
                                       const TacsScalar Xpts[],
                                       const TacsScalar vars[],
                                       const TacsScalar dvars[],
                                       const TacsScalar ddvars[],
                                       int dvLen, TacsScalar dfdx[] ){
  // Retrieve the basis object for this element (if it exists)
  TACSElementBasis *basis = element->getElementBasis();

  if (basis){
    for ( int i = 0; i < basis->getNumQuadraturePoints(); i++ ){
      double pt[3];
      double weight = basis->getQuadraturePoint(i, pt);

      // Evaluate the strain energy density
      TacsScalar U0 = 0.0;
      int count = element->evalPointQuantity(elemIndex, TACS_STRAIN_ENERGY_DENSITY,
                                             time, i, pt,
                                             Xpts, vars, dvars, ddvars, &U0);

      if (count >= 1){
        // Evaluate the determinant of the Jacobian
        TacsScalar Xd[9], J[9];
        TacsScalar detJ = basis->getJacobianTransform(pt, Xpts, Xd, J);

        // Evaluate the derivative of the strain energy
        TacsScalar dfdq = detJ*weight;
        element->addPointQuantityDVSens(elemIndex, TACS_STRAIN_ENERGY_DENSITY,
                                        time, scale, i, pt,
                                        Xpts, vars, dvars, ddvars,
                                        &dfdq, dvLen, dfdx);
      }
    }
  }
}