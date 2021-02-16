/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsMatrixSlice.h
 * @brief Provides a uniform interface to iterate rows and columns in different
 * underlying matrix storage formats
 * @author Leona Gottwald
 */

#include <cstddef>
#include <iterator>
#include <vector>

#ifndef UTIL_HIGHS_MATRIX_SLICE_H_
#define UTIL_HIGHS_MATRIX_SLICE_H_

template <typename StorageFormat>
class HighsMatrixSlice;

struct HighsEmptySlice : public HighsMatrixSlice<HighsEmptySlice> {
  using HighsMatrixSlice<HighsEmptySlice>::HighsMatrixSlice;
};
struct HighsCompressedSlice : public HighsMatrixSlice<HighsCompressedSlice> {
  using HighsMatrixSlice<HighsCompressedSlice>::HighsMatrixSlice;
};
struct HighsIndexedSlice : public HighsMatrixSlice<HighsIndexedSlice> {
  using HighsMatrixSlice<HighsIndexedSlice>::HighsMatrixSlice;
};
struct HighsTripletListSlice : public HighsMatrixSlice<HighsTripletListSlice> {
  using HighsMatrixSlice<HighsTripletListSlice>::HighsMatrixSlice;
};
struct HighsTripletTreeSliceInOrder
    : public HighsMatrixSlice<HighsTripletTreeSliceInOrder> {
  using HighsMatrixSlice<HighsTripletTreeSliceInOrder>::HighsMatrixSlice;
};
struct HighsTripletTreeSlicePreOrder
    : public HighsMatrixSlice<HighsTripletTreeSlicePreOrder> {
  using HighsMatrixSlice<HighsTripletTreeSlicePreOrder>::HighsMatrixSlice;
};
struct HighsTripletPositionSlice
    : public HighsMatrixSlice<HighsTripletPositionSlice> {
  using HighsMatrixSlice<HighsTripletPositionSlice>::HighsMatrixSlice;
};

class HighsSliceNonzero {
  template <typename>
  friend class HighsMatrixSlice;

 private:
  const int* index_;
  const double* value_;

 public:
  HighsSliceNonzero() = default;
  HighsSliceNonzero(const int* index, const double* value)
      : index_(index), value_(value) {}
  int index() const { return *index_; }
  double value() const { return *value_; }
};

template <>
class HighsMatrixSlice<HighsEmptySlice> {
  using iterator = int;

  static constexpr int begin() { return 0; }
  static constexpr int end() { return 0; }
};

template <>
class HighsMatrixSlice<HighsCompressedSlice> {
  const int* index;
  const double* value;
  int len;

 public:
  class iterator {
    HighsSliceNonzero pos_;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(const int* index, const double* value) : pos_(index, value) {}
    iterator() = default;

    iterator operator++(int) {
      iterator prev = *this;
      ++pos_.index_;
      ++pos_.value_;
      return prev;
    }
    iterator& operator++() {
      ++pos_.index_;
      ++pos_.value_;
      return *this;
    }
    reference operator*() const { return pos_; }
    pointer operator->() const { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;
      i.pos_.index_ += v;
      i.pos_.value_ += v;
      return i;
    }

    bool operator==(const iterator& rhs) const {
      return pos_.index_ == rhs.pos_.index_;
    }
    bool operator!=(const iterator& rhs) const {
      return pos_.index_ != rhs.pos_.index_;
    }
  };

  HighsMatrixSlice(const int* index, const double* value, int len)
      : index(index), value(value), len(len) {}
  iterator begin() const { return iterator{index, value}; }
  iterator end() const { return iterator{index + len, nullptr}; }
};

template <>
class HighsMatrixSlice<HighsIndexedSlice> {
  const int* index;
  const double* denseValues;
  int len;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const double* denseValues;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(const int* index, const double* denseValues)
        : pos_(index, denseValues), denseValues(denseValues) {}
    iterator() = default;

    iterator operator++(int) {
      iterator prev = *this;
      ++pos_.index_;
      return prev;
    }
    iterator& operator++() {
      ++pos_.index_;
      return *this;
    }
    reference operator*() {
      pos_.value_ = &denseValues[*pos_.index_];
      return pos_;
    }
    pointer operator->() {
      pos_.value_ = &denseValues[*pos_.index_];
      return &pos_;
    }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    bool operator==(const iterator& rhs) const {
      return pos_.index_ == rhs.pos_.index_;
    }
    bool operator!=(const iterator& rhs) const {
      return pos_.index_ != rhs.pos_.index_;
    }
  };

  HighsMatrixSlice(const int* index, const double* denseValues, int len)
      : index(index), denseValues(denseValues), len(len) {}
  iterator begin() const { return iterator{index, denseValues}; }
  iterator end() const { return iterator{index + len, nullptr}; }
};

template <>
class HighsMatrixSlice<HighsTripletListSlice> {
  const int* nodeIndex;
  const double* nodeValue;
  const int* nodeNext;
  int head;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const int* nodeNext;
    int currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(int node) : currentNode(node) {}
    iterator(const int* nodeIndex, const double* nodeValue, const int* nodeNext,
             int node)
        : pos_(nodeIndex + node, nodeValue + node),
          nodeNext(nodeNext),
          currentNode(node) {}
    iterator() = default;

    iterator& operator++() {
      pos_.index_ -= currentNode;
      pos_.value_ -= currentNode;
      currentNode = nodeNext[currentNode];
      pos_.index_ += currentNode;
      pos_.value_ += currentNode;
      return *this;
    }
    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() { return pos_; }
    pointer operator->() { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    int position() const { return currentNode; }

    bool operator==(const iterator& rhs) const {
      return currentNode == rhs.currentNode;
    }
    bool operator!=(const iterator& rhs) const {
      return currentNode != rhs.currentNode;
    }
  };

  HighsMatrixSlice(const int* nodeIndex, const double* nodeValue,
                   const int* nodeNext, int head)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodeNext(nodeNext),
        head(head) {}
  iterator begin() const {
    return iterator{nodeIndex, nodeValue, nodeNext, head};
  }
  iterator end() const { return iterator{-1}; }
};

template <>
class HighsMatrixSlice<HighsTripletTreeSlicePreOrder> {
  const int* nodeIndex;
  const double* nodeValue;
  const int* nodeLeft;
  const int* nodeRight;
  int root;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const int* nodeLeft;
    const int* nodeRight;
    std::vector<int> stack;
    int currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(int node) : currentNode(node) {}
    iterator(const int* nodeIndex, const double* nodeValue, const int* nodeLeft,
             const int* nodeRight, int node)
        : pos_(nodeIndex + node, nodeValue + node),
          nodeLeft(nodeLeft),
          nodeRight(nodeRight),
          currentNode(node) {
      stack.reserve(16);
      stack.push_back(-1);
    }
    iterator() = default;

    iterator& operator++() {
      int offset = -currentNode;
      if (nodeLeft[currentNode] != -1) {
        if (nodeRight[currentNode] != -1)
          stack.push_back(nodeRight[currentNode]);
        currentNode = nodeLeft[currentNode];
      } else {
        currentNode = stack.back();
        stack.pop_back();
      }
      offset += currentNode;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return *this;
    }

    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() { return pos_; }
    pointer operator->() { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    int position() const { return currentNode; }

    bool operator==(const iterator& rhs) const {
      return currentNode == rhs.currentNode;
    }
    bool operator!=(const iterator& rhs) const {
      return currentNode != rhs.currentNode;
    }
  };

  HighsMatrixSlice(const int* nodeIndex, const double* nodeValue,
                   const int* nodeLeft, const int* nodeRight, int root)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodeLeft(nodeLeft),
        nodeRight(nodeRight),
        root(root) {}

  iterator end() const { return iterator{-1}; }
  iterator begin() const {
    // avoid allocation if there are no elements
    if (root == -1) return end();
    return iterator{nodeIndex, nodeValue, nodeLeft, nodeRight, root};
  }
};

template <>
class HighsMatrixSlice<HighsTripletTreeSliceInOrder> {
  const int* nodeIndex;
  const double* nodeValue;
  const int* nodeLeft;
  const int* nodeRight;
  int root;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const int* nodeLeft;
    const int* nodeRight;
    std::vector<int> stack;
    int currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(int node) : currentNode(node) {}
    iterator(const int* nodeIndex, const double* nodeValue, const int* nodeLeft,
             const int* nodeRight, int node)
        : pos_(nodeIndex, nodeValue),
          nodeLeft(nodeLeft),
          nodeRight(nodeRight),
          currentNode(node) {
      stack.reserve(16);
      stack.push_back(-1);
      if (currentNode == -1) return;
      while (nodeLeft[currentNode] != -1) {
        stack.push_back(currentNode);
        currentNode = nodeLeft[currentNode];
      }

      pos_.index_ += currentNode;
      pos_.value_ += currentNode;
    }
    iterator() = default;

    iterator& operator++() {
      int offset = -currentNode;
      if (nodeRight[currentNode] != -1) {
        currentNode = nodeRight[currentNode];
        while (nodeLeft[currentNode] != -1) {
          stack.push_back(currentNode);
          currentNode = nodeLeft[currentNode];
        }
      } else {
        currentNode = stack.back();
        stack.pop_back();
      }
      offset += currentNode;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return *this;
    }

    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() { return pos_; }
    pointer operator->() { return &pos_; }
    iterator operator+(difference_type v) const {
      iterator i = *this;

      while (v > 0) {
        --v;
        ++i;
      }

      return i;
    }

    int position() const { return currentNode; }

    bool operator==(const iterator& rhs) const {
      return currentNode == rhs.currentNode;
    }
    bool operator!=(const iterator& rhs) const {
      return currentNode != rhs.currentNode;
    }
  };

  HighsMatrixSlice(const int* nodeIndex, const double* nodeValue,
                   const int* nodeLeft, const int* nodeRight, int root)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodeLeft(nodeLeft),
        nodeRight(nodeRight),
        root(root) {}

  iterator end() const { return iterator{-1}; }
  iterator begin() const {
    // avoid allocation if there are no elements
    if (root == -1) return end();
    return iterator{nodeIndex, nodeValue, nodeLeft, nodeRight, root};
  }
};

template <>
class HighsMatrixSlice<HighsTripletPositionSlice> {
  const int* nodeIndex;
  const double* nodeValue;
  const int* nodePositions;
  int len;

 public:
  class iterator {
    HighsSliceNonzero pos_;
    const int* node;
    int currentNode;

   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HighsSliceNonzero;
    using difference_type = std::ptrdiff_t;
    using pointer = const HighsSliceNonzero*;
    using reference = const HighsSliceNonzero&;

    iterator(const int* node) : node(node) {}
    iterator(const int* nodeIndex, const double* nodeValue, const int* node)
        : pos_(nodeIndex, nodeValue), currentNode(0) {}
    iterator() = default;

    iterator& operator++() {
      ++node;
      return *this;
    }

    iterator operator++(int) {
      iterator prev = *this;
      ++(*this);
      return prev;
    }
    reference operator*() {
      int offset = -currentNode + *node;
      currentNode = *node;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return pos_;
    }
    pointer operator->() {
      int offset = -currentNode + *node;
      currentNode = *node;
      pos_.index_ += offset;
      pos_.value_ += offset;
      return &pos_;
    }
    iterator operator+(difference_type v) const {
      iterator i = *this;
      i.node += v;
      return i;
    }

    int position() const { return currentNode; }

    bool operator==(const iterator& rhs) const {
      return currentNode == rhs.currentNode;
    }

    bool operator!=(const iterator& rhs) const {
      return currentNode != rhs.currentNode;
    }
  };

  HighsMatrixSlice(const int* nodeIndex, const double* nodeValue,
                   const int* nodePositions, int len)
      : nodeIndex(nodeIndex),
        nodeValue(nodeValue),
        nodePositions(nodePositions),
        len(len) {}

  iterator begin() { return iterator{nodeIndex, nodeValue, nodePositions}; }

  iterator end() const { return iterator{nodePositions + len}; }
};

#endif
