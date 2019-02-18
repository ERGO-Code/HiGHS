/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file interfaces/OsiHiGHSInterface.hpp
 * @brief Osi/HiGHS interface header
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#ifndef OsiHiGHSSolverInterface_H
#define OsiHiGHSSolverInterface_H

#include "OsiSolverInterface.hpp"
#include "Highs.h"


/** HiGHS Solver Interface
 *
 *  Instantiation of OsiSolverInterface for HiGHS
 */
class OsiHiGHSSolverInterface : virtual public OsiSolverInterface {

public:
   //---------------------------------------------------------------------------
   /**@name Solve methods */
   //@{
   /// Solve initial LP relaxation
   /// @todo implement
   virtual void initialSolve() { };

   /// Resolve an LP relaxation after problem modification
   /// @todo implement
   virtual void resolve() { };

   /// Invoke solver's built-in enumeration algorithm
   /// @todo implement
   virtual void branchAndBound() { }
   //@}

   //---------------------------------------------------------------------------
   ///@name Parameter set/get methods
   ///@todo use OsiSolverInterface default implementation or override?
   ///@{
   // Set an integer parameter
   // bool setIntParam(OsiIntParam key, int value);
   // Set an double parameter
   // bool setDblParam(OsiDblParam key, double value);
   // Set a string parameter
   // bool setStrParam(OsiStrParam key, const std::string &value);
   // Get an integer parameter
   // bool getIntParam(OsiIntParam key, int &value) const;
   // Get an double parameter
   // bool getDblParam(OsiDblParam key, double &value) const;
   // Get a string parameter
   // bool getStrParam(OsiStrParam key, std::string &value) const;
   //@}

   //---------------------------------------------------------------------------
   ///@name Methods returning info on how the solution process terminated
   ///@{
   ///@todo implement
   /// Are there a numerical difficulties?
   virtual bool isAbandoned() const { return false; }
   /// Is optimality proven?
   virtual bool isProvenOptimal() const { return false; }
   /// Is primal infeasiblity proven?
   virtual bool isProvenPrimalInfeasible() const { return false; }
   /// Is dual infeasiblity proven?
   virtual bool isProvenDualInfeasible() const { return false; }
   /// Is the given primal objective limit reached?
   virtual bool isPrimalObjectiveLimitReached() const { return false; }
   /// Is the given dual objective limit reached?
   virtual bool isDualObjectiveLimitReached() const { return false; }
   /// Iteration limit reached?
   virtual bool isIterationLimitReached() const { return false; }
   //@}

   //---------------------------------------------------------------------------
   ///@name Warm start methods
   ///@{
   ///@todo implement

   /// Get an empty warm start object
   CoinWarmStart *getEmptyWarmStart() const { return NULL; }

   /// Get warmstarting information
   virtual CoinWarmStart* getWarmStart() const { return NULL; }

   /** Set warmstarting information. Return true/false depending on whether
    *  the warmstart information was accepted or not.
    */
   virtual bool setWarmStart(const CoinWarmStart *warmstart) { return false; }
   ///@}

   //---------------------------------------------------------------------------
   ///@name Problem query methods
   ///@{
   ///@todo implement

   /// Get number of columns
   virtual int getNumCols() const { return 0; }

   /// Get number of rows
   virtual int getNumRows() const { return 0; }

   /// Get number of nonzero elements
   virtual int getNumElements() const { return 0; }

   /// Get pointer to array[getNumCols()] of column lower bounds
   virtual const double* getColLower() const { return NULL; }

   /// Get pointer to array[getNumCols()] of column upper bounds
   virtual const double* getColUpper() const { return NULL; }

   /// Get pointer to array[getNumRows()] of row constraint senses.
   virtual const char* getRowSense() const { return NULL; }

   /// Get pointer to array[getNumRows()] of rows right-hand sides
   virtual const double* getRightHandSide() const { return NULL; }

   /// Get pointer to array[getNumRows()] of row ranges.
   virtual const double* getRowRange() const { return NULL; }

   /// Get pointer to array[getNumRows()] of row lower bounds
   virtual const double* getRowLower() const { return NULL; }

   /// Get pointer to array[getNumRows()] of row upper bounds
   virtual const double* getRowUpper() const { return NULL; }

   /// Get pointer to array[getNumCols()] of objective function coefficients
   virtual const double* getObjCoefficients() const { return NULL; }

   /// Get objective function sense (1 for min (default), -1 for max)
   virtual double getObjSense() const  { return 1; }

   /// Return true if column is continuous
   virtual bool isContinuous(int colNumber) const  { return true; }

   /// Get pointer to row-wise copy of matrix
   virtual const CoinPackedMatrix* getMatrixByRow() const { return NULL; }

   /// Get pointer to column-wise copy of matrix
   virtual const CoinPackedMatrix *getMatrixByCol() const { return NULL; }

   /// Get solver's value for infinity
   virtual double getInfinity() const { return 100.0; }
   //@}

   ///@name Solution query methods
   ///@{
   ///@todo implement

   /// Get pointer to array[getNumCols()] of primal solution vector
   virtual const double* getColSolution() const { return NULL; }

   /// Get pointer to array[getNumRows()] of dual prices
   virtual const double* getRowPrice() const { return NULL; }

   /// Get a pointer to array[getNumCols()] of reduced costs
   virtual const double* getReducedCost() const { return NULL; }

   /// Get pointer to array[getNumRows()] of row activity levels (constraint matrix times the solution vector)
   virtual const double *getRowActivity() const { return NULL; }

   /// Get objective function value
   virtual double getObjValue() const { return 0.0; }

   /// Get how many iterations it took to solve the problem (whatever "iteration" mean to the solver)
   virtual int getIterationCount() const { return 42; }

   /// Get as many dual rays as the solver can provide.
   virtual std::vector< double*> getDualRays(int maxNumRays, bool fullRay = false) const { return std::vector<double*>(0); }

   /// Get as many primal rays as the solver can provide.
   virtual std::vector<double*> getPrimalRays(int maxNumRays) const { return std::vector<double*>(0); }

   ///@}

//   //---------------------------------------------------------------------------

//   /**@name Problem modifying methods */
//   //@{
//   //-------------------------------------------------------------------------
//   /**@name Changing bounds on variables and constraints */
//   //@{
//   /** Set an objective function coefficient */
//   virtual void setObjCoeff(int elementIndex, double elementValue);

//   /** Set a a set of objective function coefficients */
//   virtual void setObjCoeffSet(const int *indexFirst,
//     const int *indexLast,
//     const double *coeffList);

//   using OsiSolverInterface::setColLower;
//   /** Set a single column lower bound<br>
//     	  Use -COIN_DBL_MAX for -infinity. */
//   virtual void setColLower(int elementIndex, double elementValue);

//   using OsiSolverInterface::setColUpper;
//   /** Set a single column upper bound<br>
//     	  Use COIN_DBL_MAX for infinity. */
//   virtual void setColUpper(int elementIndex, double elementValue);

//   /** Set a single column lower and upper bound<br>
//     	  The default implementation just invokes <code>setColLower()</code> and
//     	  <code>setColUpper()</code> */
//   virtual void setColBounds(int elementIndex,
//     double lower, double upper);

//   /** Set the bounds on a number of columns simultaneously<br>
//     	  The default implementation just invokes <code>setCollower()</code> and
//     	  <code>setColupper()</code> over and over again.
//     	  @param <code>[indexfirst,indexLast]</code> contains the indices of
//     	         the constraints whose </em>either</em> bound changes
//     	  @param boundList the new lower/upper bound pairs for the variables
//       */
//   virtual void setColSetBounds(const int *indexFirst,
//     const int *indexLast,
//     const double *boundList);

//   /** Set a single row lower bound<br>
//     	  Use -COIN_DBL_MAX for -infinity. */
//   virtual void setRowLower(int elementIndex, double elementValue);

//   /** Set a single row upper bound<br>
//     	  Use COIN_DBL_MAX for infinity. */
//   virtual void setRowUpper(int elementIndex, double elementValue);

//   /** Set a single row lower and upper bound<br>
//     	  The default implementation just invokes <code>setRowLower()</code> and
//     	  <code>setRowUpper()</code> */
//   virtual void setRowBounds(int elementIndex,
//     double lower, double upper);

//   /** Set the type of a single row<br> */
//   virtual void setRowType(int index, char sense, double rightHandSide,
//     double range);

//   /** Set the bounds on a number of rows simultaneously<br>
//     	  The default implementation just invokes <code>setRowLower()</code> and
//     	  <code>setRowUpper()</code> over and over again.
//     	  @param <code>[indexfirst,indexLast]</code> contains the indices of
//     	         the constraints whose </em>either</em> bound changes
//     	  @param boundList the new lower/upper bound pairs for the constraints
//       */
//   virtual void setRowSetBounds(const int *indexFirst,
//     const int *indexLast,
//     const double *boundList);

//   /** Set the type of a number of rows simultaneously<br>
//     	  The default implementation just invokes <code>setRowType()</code> and
//     	  over and over again.
//     	  @param <code>[indexfirst,indexLast]</code> contains the indices of
//     	         the constraints whose type changes
//     	  @param senseList the new senses
//     	  @param rhsList   the new right hand sides
//     	  @param rangeList the new ranges
//       */
//   virtual void setRowSetTypes(const int *indexFirst,
//     const int *indexLast,
//     const char *senseList,
//     const double *rhsList,
//     const double *rangeList);
//   //@}

//   //-------------------------------------------------------------------------
//   /**@name Integrality related changing methods */
//   //@{
//   /** Set the index-th variable to be a continuous variable */
//   virtual void setContinuous(int index);
//   /** Set the index-th variable to be an integer variable */
//   virtual void setInteger(int index);
//   /** Set the variables listed in indices (which is of length len) to be
// 	  continuous variables */
//   virtual void setContinuous(const int *indices, int len);
//   /** Set the variables listed in indices (which is of length len) to be
// 	  integer variables */
//   virtual void setInteger(const int *indices, int len);
//   //@}

//   //-------------------------------------------------------------------------
//   /// Set objective function sense (1 for min (default), -1 for max,)
//   virtual void setObjSense(double s);

//   /** Set the primal solution column values
    
//     	colsol[numcols()] is an array of values of the problem column
//     	variables. These values are copied to memory owned by the
//     	solver object or the solver.  They will be returned as the
//     	result of colsol() until changed by another call to
//     	setColsol() or by a call to any solver routine.  Whether the
//     	solver makes use of the solution in any way is
//     	solver-dependent. 
//     */
//   virtual void setColSolution(const double *colsol);

//   /** Set dual solution vector
    
//     	rowprice[numrows()] is an array of values of the problem row
//     	dual variables. These values are copied to memory owned by the
//     	solver object or the solver.  They will be returned as the
//     	result of rowprice() until changed by another call to
//     	setRowprice() or by a call to any solver routine.  Whether the
//     	solver makes use of the solution in any way is
//     	solver-dependent. 
//     */
//   virtual void setRowPrice(const double *rowprice);

//   //-------------------------------------------------------------------------
//   /**@name Methods to expand a problem.<br>
//        Note that if a column is added then by default it will correspond to a
//        continuous variable. */
//   //@{
//   using OsiSolverInterface::addCol;
//   /** */
//   virtual void addCol(const CoinPackedVectorBase &vec,
//     const double collb, const double colub,
//     const double obj);

//   using OsiSolverInterface::addCols;
//   /** */
//   virtual void addCols(const int numcols,
//     const CoinPackedVectorBase *const *cols,
//     const double *collb, const double *colub,
//     const double *obj);
//   /** */
//   virtual void deleteCols(const int num, const int *colIndices);

//   using OsiSolverInterface::addRow;
//   /** */
//   virtual void addRow(const CoinPackedVectorBase &vec,
//     const double rowlb, const double rowub);
//   /** */
//   virtual void addRow(const CoinPackedVectorBase &vec,
//     const char rowsen, const double rowrhs,
//     const double rowrng);

//   using OsiSolverInterface::addRows;
//   /** */
//   virtual void addRows(const int numrows,
//     const CoinPackedVectorBase *const *rows,
//     const double *rowlb, const double *rowub);
//   /** */
//   virtual void addRows(const int numrows,
//     const CoinPackedVectorBase *const *rows,
//     const char *rowsen, const double *rowrhs,
//     const double *rowrng);
//   /** */
//   virtual void deleteRows(const int num, const int *rowIndices);
//   //@}
//   //@}

//   //---------------------------------------------------------------------------

//   /**@name Methods to input a problem */
//   //@{
//   /** Load in an problem by copying the arguments (the constraints on the
//         rows are given by lower and upper bounds). If a pointer is 0 then the
//         following values are the default:
//         <ul>
//           <li> <code>colub</code>: all columns have upper bound infinity
//           <li> <code>collb</code>: all columns have lower bound 0 
//           <li> <code>rowub</code>: all rows have upper bound infinity
//           <li> <code>rowlb</code>: all rows have lower bound -infinity
// 	  <li> <code>obj</code>: all variables have 0 objective coefficient
//         </ul>
//     */
//   virtual void loadProblem(const CoinPackedMatrix &matrix,
//     const double *collb, const double *colub,
//     const double *obj,
//     const double *rowlb, const double *rowub);

//   /** Load in an problem by assuming ownership of the arguments (the
//         constraints on the rows are given by lower and upper bounds). For
//         default values see the previous method. <br>
// 	<strong>WARNING</strong>: The arguments passed to this method will be
// 	freed using the C++ <code>delete</code> and <code>delete[]</code>
// 	functions. 
//     */
//   virtual void assignProblem(CoinPackedMatrix *&matrix,
//     double *&collb, double *&colub, double *&obj,
//     double *&rowlb, double *&rowub);

//   /** Load in an problem by copying the arguments (the constraints on the
// 	rows are given by sense/rhs/range triplets). If a pointer is 0 then the
// 	following values are the default:
// 	<ul>
//           <li> <code>colub</code>: all columns have upper bound infinity
//           <li> <code>collb</code>: all columns have lower bound 0 
// 	  <li> <code>obj</code>: all variables have 0 objective coefficient
//           <li> <code>rowsen</code>: all rows are >=
//           <li> <code>rowrhs</code>: all right hand sides are 0
//           <li> <code>rowrng</code>: 0 for the ranged rows
//         </ul>
//     */
//   virtual void loadProblem(const CoinPackedMatrix &matrix,
//     const double *collb, const double *colub,
//     const double *obj,
//     const char *rowsen, const double *rowrhs,
//     const double *rowrng);

//   /** Load in an problem by assuming ownership of the arguments (the
//         constraints on the rows are given by sense/rhs/range triplets). For
//         default values see the previous method. <br>
// 	<strong>WARNING</strong>: The arguments passed to this method will be
// 	freed using the C++ <code>delete</code> and <code>delete[]</code>
// 	functions. 
//     */
//   virtual void assignProblem(CoinPackedMatrix *&matrix,
//     double *&collb, double *&colub, double *&obj,
//     char *&rowsen, double *&rowrhs,
//     double *&rowrng);

//   /** Just like the other loadProblem() methods except that the matrix is
// 	given in a standard column major ordered format (without gaps). */
//   virtual void loadProblem(const int numcols, const int numrows,
//     const int *start, const int *index,
//     const double *value,
//     const double *collb, const double *colub,
//     const double *obj,
//     const double *rowlb, const double *rowub);

//   /** Just like the other loadProblem() methods except that the matrix is
// 	given in a standard column major ordered format (without gaps). */
//   virtual void loadProblem(const int numcols, const int numrows,
//     const int *start, const int *index,
//     const double *value,
//     const double *collb, const double *colub,
//     const double *obj,
//     const char *rowsen, const double *rowrhs,
//     const double *rowrng);

//   using OsiSolverInterface::readMps;
//   /** Read an mps file from the given filename */
//   virtual int readMps(const char *filename,
//     const char *extension = "mps");

//   /** Write the problem into an mps file of the given filename.
// 	If objSense is non zero then -1.0 forces the code to write a
// 	maximization objective and +1.0 to write a minimization one.
// 	If 0.0 then solver can do what it wants */
//   virtual void writeMps(const char *filename,
//     const char *extension = "mps",
//     double objSense = 0.0) const;

//   //@}

//   /**@name Message handling */
//   //@{
//   /** Pass in a message handler
//         It is the client's responsibility to destroy a message handler installed
//         by this routine; it will not be destroyed when the solver interface is
//         destroyed. 
//      */
//   void passInMessageHandler(CoinMessageHandler *handler);
//   //@}

//   /**@name Constructors and destructor */
//   //@{
//   /// Default Constructor
  OsiHiGHSSolverInterface();

//   /// Clone
//   virtual OsiSolverInterface *clone(bool copyData = true) const;

//   /// Copy constructor
//   OsiHiGHSSolverInterface(const OsiHiGHSSolverInterface &);

//   /// Assignment operator
//   OsiHiGHSSolverInterface &operator=(const OsiHiGHSSolverInterface &rhs);

  /// Destructor
  virtual ~OsiHiGHSSolverInterface();

//   /// Resets as if default constructor
//   virtual void reset();
//   //@}

//   /***************************************************************************/
//   /**@name OsiSimplexInterface methods 
//   */
//   //@{

//   /** Returns 1 if can just do getBInv etc
//       2 if has all OsiSimplex methods
//       and 0 if it has none */
//   virtual int canDoSimplexInterface() const;

//   using OsiSolverInterface::enableSimplexInterface;
//   /** Useless function, defined only for compatibility with 
//      OsiSimplexInterface
//   */
//   virtual void enableSimplexInterface(int doingPrimal) {};

//   /** Useless function, defined only for compatibility with 
//      OsiSimplexInterface
//   */
//   virtual void disableSimplexInterface() {};

//   /** Useless function, defined only for compatibility with 
//      OsiSimplexInterface
//   */
//   virtual void enableFactorization() const {};

//   /** Useless function, defined only for compatibility with 
//      OsiSimplexInterface
//   */
//   virtual void disableFactorization() const {};

//   ///Returns true if a basis is available
//   virtual bool basisIsAvailable() const;

//   /** Returns a basis status of the structural/artificial variables 
//      At present as warm start i.e 0: free, 1: basic, 2: upper, 3: lower
//   */
//   virtual void getBasisStatus(int *cstat, int *rstat) const;

//   ///Get a row of the tableau (slack part in slack if not NULL)
//   virtual void getBInvARow(int row, double *z, double *slack = NULL) const;

//   ///Get a row of the basis inverse
//   virtual void getBInvRow(int row, double *z) const;

//   ///Get a column of the tableau
//   virtual void getBInvACol(int col, double *vec) const;

//   ///Get a column of the basis inverse
//   virtual void getBInvCol(int col, double *vec) const;

//   /**  Get indices of the pivot variable in each row
//       (order of indices corresponds to the
//       order of elements in a vector retured by getBInvACol() and
//       getBInvCol()).
//   */
//   virtual void getBasics(int *index) const;

//   //@}
//   /***************************************************************************/

// protected:
//   /**@name Protected methods */
//   //@{
//   /// Apply a row cut. Return true if cut was applied.
//   virtual void applyRowCut(const OsiRowCut &rc);

//   /** Apply a column cut (bound adjustment). 
//       Return true if cut was applied.
//   */
//   virtual void applyColCut(const OsiColCut &cc);
//   //@}


private:

  Highs* highs;

};

#endif
