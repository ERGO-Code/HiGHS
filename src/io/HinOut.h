/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/HinOut.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IO_HINOUT_H_
#define IO_HINOUT_H_

#include <cstring>
#include "HModel.h"
#include "Presolve.h"

class HinOut : public HPreData {
 public:
  int AcountX;
  string fileIn, fileOut;

  HinOut(string filenameIn, string filenameOut);
  void getData(HModel& ptr_model);
  void setData(HModel& ptr_model);
  void readDataPostsolve(HModel& ptr_model);
  void readDataColumnWise();
  void writeDataColumnWise();
  void clearData();
  void HinOutTestIO(HModel& ptr);
  void HinOutTestRead(HModel& ptr);
  void HinOutTestWrite(HModel& ptr);
  void compareData(int lvl);
  double getdiff(double v1, double v2);

  // data we are getting and printing
  int onumCol;
  int onumRow;

  vector<int> oAstart;
  vector<int> oAindex;
  vector<double> oAvalue;

  vector<double> ocolCost;
  vector<double> ocolLower;
  vector<double> ocolUpper;
  vector<double> orowLower;
  vector<double> orowUpper;

  int oAcountX;
};

#endif /* IO_HINOUT_H_ */
