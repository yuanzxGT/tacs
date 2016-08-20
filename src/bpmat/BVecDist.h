#ifndef TACS_BVEC_DISTRIBUTE_H
#define TACS_BVEC_DISTRIBUTE_H

/*
  Distribute/gather elements of a BVec

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "TACSObject.h"

enum TACSBVecOperation { INSERT_VALUES, ADD_VALUES };

/*
  Declare the TACSBVecDistCtx class
*/
class TACSBVecDistCtx;

/*!
  Variable map for the parallel distribution of a vector

  This class defines the mapping between the variables and processors
  and should be instantiated once for each analysis model.
*/
class TACSVarMap : public TACSObject {
 public:
  TACSVarMap( MPI_Comm _comm, int _N );
  ~TACSVarMap();

  int getDim();
  MPI_Comm getMPIComm();
  void getOwnerRange( const int **_ownerRange );

 private:
  MPI_Comm comm; // The MPI communicator
  int *ownerRange; // The ownership range of the variables
  int N; // Number of nodes on this processor
};

/*
  A class containing a pointer to an array of indices.  These indices
  may be either sorted in asscending order and unique or completely
  arbitrary. These cases are detected automatically by the index
  class. The indices may not be negative.  Negative indices are
  replaced by 0, so be warned!!

  This object takes ownership of the array that is passed in.
*/
class TACSBVecIndices : public TACSObject {
 public:
  TACSBVecIndices( int **_indices, int _nindices );
  ~TACSBVecIndices();

  // Retrieve information about the indices
  // --------------------------------------
  int getNumIndices();
  int getIndices( const int **_indices );
  int isSorted();

  // Set up/use an arg-sorted array to find the reverse
  // map to find k such that indices[k] = var 
  // --------------------------------------------------
  void setUpInverse();
  int findIndex( int index );

 private:
  int *indices;
  int nindices;
  int issorted;
  int *index_args;
};

/*!
  Distribute vector components to other processors and collect 
  contributions from other processors.

  This class is used to pass external interface variables between
  processes during parallel matrix-vector products.

  This class performs the following operation:

  for i = 1:nvars:
  local[i] = vec[vars[i]]

  where vars[i] are possibly non-local variables.

  Additionally, a reverse communication is also permitted where the
  following operation takes place,

  for i = 1:nvars:
  vec[vars[i]] += local[i]

  This operation is useful for assembling the residual equations
  within the finite--element method.
*/
class TACSBVecDistribute : public TACSObject {
 public:
  TACSBVecDistribute( TACSVarMap *rmap, TACSBVecIndices *bindex );
  ~TACSBVecDistribute();

  // Create a context to send/recv the data
  // --------------------------------------
  TACSBVecDistCtx *createCtx( int bsize );

  // Get the size of the local array
  // All arrays passed must be at least this size
  // --------------------------------------------
  int getDim();
  TACSBVecIndices *getIndices();

  // Transfer the data to the array provided
  // ---------------------------------------
  void beginForward( TACSBVecDistCtx *ctx,
                     TacsScalar *global, TacsScalar *local );
  void endForward( TACSBVecDistCtx *ctx,
                   TacsScalar *global, TacsScalar *local );

  // Add or insert data back into the vector
  // ---------------------------------------
  void beginReverse( TACSBVecDistCtx *ctx,
                     TacsScalar *local, TacsScalar *global, 
                     TACSBVecOperation op=ADD_VALUES );
  void endReverse( TACSBVecDistCtx *ctx,
                   TacsScalar *local, TacsScalar *global,
                   TACSBVecOperation op=ADD_VALUES );

  MPI_Comm getMPIComm();
  const char *TACSObjectName();

 private:
  // Block-specific implementation pointers
  // --------------------------------------
  void initImpl( int bsize );
  void (*bgetvars)( int bsize, int nvars, const int *vars, int lower,
		    TacsScalar *x, TacsScalar *y, 
		    TACSBVecOperation op );
  void (*bsetvars)( int bsize, int nvars, const int *vars, int lower,
		    TacsScalar *x, TacsScalar *y, 
		    TACSBVecOperation op );

  // The communicator and the MPI data
  MPI_Comm comm;

  // Data defining the distribution of the variables
  TACSVarMap *rmap;

  // Object containing the indices of the external variables
  // -------------------------------------------------------
  TACSBVecIndices *bindex;

  // Data for dealing with indices that are not sorted or unique
  // -----------------------------------------------------------
  int sorted_flag;
  int nvars_unsorted;
  int *ext_sorted;
  int *ext_unsorted;

  // Data for collecting external variables
  // --------------------------------------
  int n_ext_proc; // number of external procs
  int *ext_proc; //  Externall processes
  int *ext_ptr; // Displacements into the local external array
  int *ext_count; // Count per proc
  int next_vars; // Number of external vars
  const int *ext_vars; // External variables requested by this process

  // Data for the requested values
  int n_req_proc; // Processes to send non-zero mesages to
  int *req_proc; // The processors to send
  int *req_ptr;  // Displacement into requested array
  int *req_count; // number of nodes to send to each proc
  int *req_vars; // Variables that have been requested

  // The size of the receiving data on this processor
  int ext_self_ptr;
  int ext_self_count;

  // Set the name of the object
  static const char *name;
};

/*
  Context for the distribute vector class. Note that this class can
  only be created by the TACSBVecDistribute class itself so that the
  data is allocated correctly.
*/
class TACSBVecDistCtx : public TACSObject {
 public:
  ~TACSBVecDistCtx();

 private:
  TACSBVecDistCtx( TACSBVecDistribute *_me, 
                   int _bsize );

  // The block size for this context
  int bsize;

  // Pointer to ensure the context is used correctly
  TACSBVecDistribute *me;
  
  // The external data sorted
  TacsScalar *ext_sorted_vals;

  // The requested values
  TacsScalar *reqvals;

  // The MPI requests for either the sends or recvs
  MPI_Request *sends;
  MPI_Request *recvs;

  // Set the send and recv tags
  int ctx_tag;
  static int tag_value;

  // Friend class declaration
  friend class TACSBVecDistribute;
};

#endif
