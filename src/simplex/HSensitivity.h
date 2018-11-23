/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSensitivity.h
 * @brief Sensitivity and ranging information for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HSENSITIVITY_H_
#define SIMPLEX_HSENSITIVITY_H_

#include <vector>
using namespace std;

// LP model size
int numCol;
int numRow;
int numTotal;

// SOBJ data
vector<double> c_up_c;
vector<double> c_up_f;
vector<int> c_up_e;
vector<int> c_up_l;
vector<double> c_dn_c;
vector<double> c_dn_f;
vector<int> c_dn_e;
vector<int> c_dn_l;

// SBND data
vector<double> b_up_b;
vector<double> b_up_f;
vector<int> b_up_e;
vector<int> b_up_l;
vector<double> b_dn_b;
vector<double> b_dn_f;
vector<int> b_dn_e;
vector<int> b_dn_l;

/**
 * @brief Sensitivity and ranging information for HiGHS  
 */
class HSensitivity {
 public:
  /**
   * @brief Compute sensitivity and ranging information
   */
  int getSensitivityData(
		      HModel *model //!< Instance of HModel class for which sensitivity and ranging data are to be generated
		      );

  /**
   * @brief Check sensitivity and ranging information
   */
  int checkSensitivityData(
		      HModel *model //!< Instance of HModel class for which sensitivity and ranging data are to be generated
		      );
 private:
  void checkSensitivityDataSolve(
				 HModel *model //!< Instance of HModel class for which sensitivity and ranging data are to be generated
				 );
};

#endif /* SIMPLEX_HSENSITIVITY_H_ */
