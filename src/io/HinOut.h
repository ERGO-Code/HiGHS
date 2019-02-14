/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
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

#include "presolve/HPreData.h"

class HinOut : public HPreData {
 public:
  int AcountX;
  string fileIn, fileOut;

  HinOut(string filenameIn, string filenameOut);
  void getData(HighsModelObject& highs_model_object);
  void setData(HighsModelObject& highs_model_object);
  void readDataPostsolve(HighsModelObject& highs_model_object);
  void readDataColumnWise();
  void writeDataColumnWise();
  void clearData();
  void HinOutTestIO(HighsModelObject& highs_model_object);
  void HinOutTestRead(HighsModelObject& highs_model_object);
  void HinOutTestWrite(HighsModelObject& highs_model_object);
  void compareData(int lvl);
  double getdiff(double v1, double v2);

  // data we are getting and printing
  int onumCol;
  int onumRow;

  std::vector<int> oAstart;
  std::vector<int> oAindex;
  std::vector<double> oAvalue;

  std::vector<double> ocolCost;
  std::vector<double> ocolLower;
  std::vector<double> ocolUpper;
  std::vector<double> orowLower;
  std::vector<double> orowUpper;

  int oAcountX;
};

#endif /* IO_HINOUT_H_ */
