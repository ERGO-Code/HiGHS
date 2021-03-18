/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsDataStack.h
 * @brief A stack of unstructured data stored as bytes
 * @author Leona Gottwald
 */

#ifndef UTIL_HIGHS_DATA_STACK_H_
#define UTIL_HIGHS_DATA_STACK_H_

#include <cstring>
#include <type_traits>
#include <vector>

class HighsDataStack {
  std::vector<char> data;
  int position;

 public:
  void resetPosition() { position = data.size(); }

  template <typename T,
            typename std::enable_if<std::is_trivially_copyable<T>::value,
                                    int>::type = 0>
  void push(const T& r) {
    int dataSize = data.size();
    data.resize(dataSize + sizeof(T));
    std::memcpy(data.data() + dataSize, &r, sizeof(T));
  }

  template <typename T,
            typename std::enable_if<std::is_trivially_copyable<T>::value,
                                    int>::type = 0>
  void pop(T& r) {
    position -= sizeof(T);
    std::memcpy(&r, data.data() + position, sizeof(T));
  }

  template <typename T>
  void push(const std::vector<T>& r) {
    int offset = data.size();
    int numData = r.size();
    // store the data
    data.resize(offset + numData * sizeof(T) + sizeof(int));
    std::memcpy(data.data() + offset, r.data(), numData * sizeof(T));
    // store the vector size
    offset += numData * sizeof(T);
    std::memcpy(data.data() + offset, &numData, sizeof(int));
  }

  template <typename T>
  void pop(std::vector<T>& r) {
    // pop the vector size
    position -= sizeof(int);
    int numData;
    std::memcpy(&numData, &data[position], sizeof(int));
    // pop the data
    position -= numData * sizeof(T);
    r.resize(numData);
    std::memcpy(r.data(), data.data() + position, numData * sizeof(T));
  }
};

#endif