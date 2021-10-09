/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_RBTREE_H_
#define HIGHS_RBTREE_H_

#include <algorithm>
#include <cassert>
#include <type_traits>

#include "util/HighsInt.h"

namespace highs {

struct RbTreeLinks {
  enum Direction {
    kLeft = 0,
    kRight = 1,
  };

  HighsInt child[2];
#ifdef HIGHSINT64
  HighsUInt parent : 63;
  HighsUInt color : 1;
#else
  HighsUInt parent : 31;
  HighsUInt color : 1;
#endif

  bool isBlack() const { return color == 0; }

  bool isRed() const { return color == 1; }

  void makeRed() { color = 1; }

  void makeBlack() { color = 0; }

  HighsUInt getColor() const { return color; }

  void setColor(HighsUInt color) { this->color = color; }

  void setParent(HighsInt p) { parent = p + 1; }

  HighsInt getParent() const { return HighsInt(parent) - 1; }
};

template <typename Impl>
struct RbTreeTraits;

template <typename Impl>
class RbTree {
  HighsInt& rootNode;
  // static_assert(
  //     std::is_same<RbTreeLinks&, decltype(static_cast<Impl*>(nullptr)
  //                                             ->getRbTreeLinks(0))>::value,
  //     "RbTree implementation must provide a getRbTreeLinks() function that "
  //     "returns a non-const reference of RbTreeLinks given the index of a
  //     node");

  // static_assert(std::is_same<const RbTreeLinks&,
  //                           decltype(static_cast<const Impl*>(nullptr)
  //                                        ->getRbTreeLinks(0))>::value,
  //              "RbTree implementation must provide a getRbTreeLinks() const "
  //              "function that returns a const reference of RbTreeLinks given
  //              " "the index of a node");
  using KeyType = typename RbTreeTraits<Impl>::KeyType;
  using Dir = RbTreeLinks::Direction;
  // static_assert(
  //    std::is_same<bool, decltype((*static_cast<KeyType*>(nullptr)) <
  //                                (*static_cast<KeyType*>(nullptr)))>::value,
  //    "RbTree implementation must provide a getKey() function that, given the
  //    " "index of a node, returns its key which must have a type that has "
  //    "operator< defined");

  bool isRed(HighsInt node) const {
    return node != -1 &&
           static_cast<const Impl*>(this)->getRbTreeLinks(node).isRed();
  }

  bool isBlack(HighsInt node) const {
    return node == -1 ||
           static_cast<const Impl*>(this)->getRbTreeLinks(node).isBlack();
  }

  void makeRed(HighsInt node) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).makeRed();
  }

  void makeBlack(HighsInt node) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).makeBlack();
  }

  void setColor(HighsInt node, HighsUInt color) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).setColor(color);
  }

  HighsUInt getColor(HighsInt node) const {
    return static_cast<const Impl*>(this)->getRbTreeLinks(node).getColor();
  }

  void setParent(HighsInt node, HighsInt parent) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).setParent(parent);
  }

  HighsInt getParent(HighsInt node) const {
    return static_cast<const Impl*>(this)->getRbTreeLinks(node).getParent();
  }

  KeyType getKey(HighsInt node) const {
    return static_cast<const Impl*>(this)->getKey(node);
  }

  static constexpr Dir opposite(Dir dir) { return Dir(1 - dir); }

  HighsInt getChild(HighsInt node, Dir dir) const {
    return static_cast<const Impl*>(this)->getRbTreeLinks(node).child[dir];
  }

  void setChild(HighsInt node, Dir dir, HighsInt child) {
    static_cast<Impl*>(this)->getRbTreeLinks(node).child[dir] = child;
  }

  void rotate(HighsInt x, Dir dir) {
    HighsInt y = getChild(x, opposite(dir));
    HighsInt yDir = getChild(y, dir);
    setChild(x, opposite(dir), yDir);
    if (yDir != -1) setParent(yDir, x);

    HighsInt pX = getParent(x);
    setParent(y, pX);

    if (pX == -1)
      rootNode = y;
    else
      setChild(pX, Dir((x != getChild(pX, dir)) ^ dir), y);

    setChild(y, dir, x);
    setParent(x, y);
  }

  void insertFixup(HighsInt z) {
    HighsInt pZ = getParent(z);
    while (isRed(pZ)) {
      HighsInt zGrandParent = getParent(pZ);
      assert(zGrandParent != -1);

      Dir dir = Dir(getChild(zGrandParent, RbTreeLinks::kLeft) == pZ);

      HighsInt y = getChild(zGrandParent, dir);
      if (isRed(y)) {
        makeBlack(pZ);
        makeBlack(y);
        makeRed(zGrandParent);
        z = zGrandParent;
      } else {
        if (z == getChild(pZ, dir)) {
          z = pZ;
          rotate(z, opposite(dir));
          pZ = getParent(z);
          zGrandParent = getParent(pZ);
          assert(zGrandParent != -1);
        }

        makeBlack(pZ);
        makeRed(zGrandParent);
        rotate(zGrandParent, dir);
      }

      pZ = getParent(z);
    }

    makeBlack(rootNode);
  }

  void transplant(HighsInt u, HighsInt v, HighsInt& nilParent) {
    HighsInt p = getParent(u);

    if (p == -1)
      rootNode = v;
    else
      setChild(p, Dir(u != getChild(p, RbTreeLinks::kLeft)), v);

    if (v == -1)
      nilParent = p;
    else
      setParent(v, p);
  }

  void deleteFixup(HighsInt x, const HighsInt nilParent) {
    while (x != rootNode && isBlack(x)) {
      Dir dir;

      HighsInt p = x == -1 ? nilParent : getParent(x);
      dir = Dir(x == getChild(p, RbTreeLinks::kLeft));
      HighsInt w = getChild(p, dir);
      assert(w != -1);

      if (isRed(w)) {
        makeBlack(w);
        makeRed(p);
        rotate(p, opposite(dir));
        assert((x == -1 && p == nilParent) || (x != -1 && p == getParent(x)));
        w = getChild(p, dir);
        assert(w != -1);
      }

      if (isBlack(getChild(w, RbTreeLinks::kLeft)) &&
          isBlack(getChild(w, RbTreeLinks::kRight))) {
        makeRed(w);
        x = p;
      } else {
        if (isBlack(getChild(w, dir))) {
          makeBlack(getChild(w, opposite(dir)));
          makeRed(w);
          rotate(w, dir);
          assert((x == -1 && p == nilParent) || (x != -1 && p == getParent(x)));
          w = getChild(p, dir);
        }
        setColor(w, getColor(p));
        makeBlack(p);
        makeBlack(getChild(w, dir));
        rotate(p, opposite(dir));
        x = rootNode;
      }
    }

    if (x != -1) makeBlack(x);
  }

 public:
  RbTree(HighsInt& rootNode) : rootNode(rootNode) {}

  bool empty() const { return rootNode == -1; }

  HighsInt first(HighsInt x) const {
    if (x == -1) return -1;

    while (true) {
      HighsInt lX = getChild(x, RbTreeLinks::kLeft);
      if (lX == -1) return x;
      x = lX;
    }
  }

  HighsInt last(HighsInt x) const {
    if (x == -1) return -1;

    while (true) {
      HighsInt rX = getChild(x, RbTreeLinks::kRight);
      if (rX == -1) return x;
      x = rX;
    }
  }

  HighsInt first() const { return first(rootNode); }

  HighsInt last() const { return last(rootNode); }

  HighsInt successor(HighsInt x) const {
    HighsInt y = getChild(x, RbTreeLinks::kRight);
    if (y != -1) return first(y);

    y = getParent(x);
    while (y != -1 && x == getChild(y, RbTreeLinks::kRight)) {
      x = y;
      y = getParent(x);
    }

    return y;
  }

  HighsInt predecessor(HighsInt x) const {
    HighsInt y = getChild(x, RbTreeLinks::kLeft);
    if (y != -1) return last(y);

    y = getParent(x);
    while (y != -1 && x == getChild(y, RbTreeLinks::kLeft)) {
      x = y;
      y = getParent(x);
    }

    return y;
  }

  std::pair<HighsInt, bool> find(const KeyType& key, HighsInt treeRoot) {
    HighsInt y = -1;
    HighsInt x = treeRoot;
    while (x != -1) {
      HighsInt cmp = 1 - (getKey(x) < key) + (key < getKey(x));
      switch (cmp) {
        case 0:
          y = x;
          x = getChild(y, RbTreeLinks::kRight);
          break;
        case 1:
          return std::make_pair(x, true);
        case 2:
          y = x;
          x = getChild(y, RbTreeLinks::kLeft);
      }
    }

    return std::make_pair(y, false);
  }

  std::pair<HighsInt, bool> find(const KeyType& key) {
    return find(key, rootNode);
  }

  void link(HighsInt z, HighsInt parent) {
    setParent(z, parent);
    if (parent == -1)
      rootNode = z;
    else
      setChild(parent, Dir(getKey(parent) < getKey(z)), z);

    setChild(z, RbTreeLinks::kLeft, -1);
    setChild(z, RbTreeLinks::kRight, -1);
    makeRed(z);
    insertFixup(z);
  }

  void link(HighsInt z) {
    HighsInt y = -1;
    HighsInt x = rootNode;
    while (x != -1) {
      y = x;
      x = getChild(y, Dir(getKey(x) < getKey(z)));
    }

    link(z, y);
  }

  void unlink(HighsInt z) {
    HighsInt nilParent = -1;
    HighsInt y = z;
    bool yWasBlack = isBlack(y);
    HighsInt x;

    if (getChild(z, RbTreeLinks::kLeft) == -1) {
      x = getChild(z, RbTreeLinks::kRight);
      transplant(z, x, nilParent);
    } else if (getChild(z, RbTreeLinks::kRight) == -1) {
      x = getChild(z, RbTreeLinks::kLeft);
      transplant(z, x, nilParent);
    } else {
      y = first(getChild(z, RbTreeLinks::kRight));
      yWasBlack = isBlack(y);
      x = getChild(y, RbTreeLinks::kRight);
      if (getParent(y) == z) {
        if (x == -1)
          nilParent = y;
        else
          setParent(x, y);
      } else {
        transplant(y, getChild(y, RbTreeLinks::kRight), nilParent);
        HighsInt zRight = getChild(z, RbTreeLinks::kRight);
        setChild(y, RbTreeLinks::kRight, zRight);
        setParent(zRight, y);
      }
      transplant(z, y, nilParent);
      HighsInt zLeft = getChild(z, RbTreeLinks::kLeft);
      setChild(y, RbTreeLinks::kLeft, zLeft);
      setParent(zLeft, y);
      setColor(y, getColor(z));
    }

    if (yWasBlack) deleteFixup(x, nilParent);
  }
};

}  // namespace highs

#endif
