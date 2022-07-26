/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_UTIL_HASH_TREE_H_
#define HIGHS_UTIL_HASH_TREE_H_

#include "util/HighsHash.h"

template <typename K, typename V = void>
class HighsHashTree {
  using Entry = HighsHashTableEntry<K, V>;
  using ValueType = typename std::remove_reference<
      decltype(reinterpret_cast<Entry*>(0x1)->value())>::type;

  enum Type {
    kEmpty = 0,
    kSingleLeaf = 1,
    kMultiLeaf = 2,
    kForwardingNode = 3,
    kBranchNode = 4,
  };

  struct SingleLeaf {
    uint64_t hash;
    HighsHashTableEntry<K, V> entry;

    template <typename... Args>
    SingleLeaf(uint64_t hash, Args&&... args)
        : hash(hash), entry(std::forward<Args>(args)...) {}
  };

  struct ListNode {
    ListNode* next;
    HighsHashTableEntry<K, V> entry;

    ListNode(HighsHashTableEntry<K, V>&& entry)
        : next(nullptr), entry(std::move(entry)) {}
  };

  struct MultiLeaf {
    uint64_t hash;
    ListNode head;

    template <typename... Args>
    MultiLeaf(uint64_t hash, Args&&... args)
        : hash(hash),
          head(HighsHashTableEntry<K, V>(std::forward<Args>(args)...)) {}
  };

  struct BranchNode;

  struct ForwardingNode;

  struct NodePtr {
    uintptr_t ptrAndType;

    NodePtr() : ptrAndType(kEmpty) {}
    NodePtr(std::nullptr_t) : ptrAndType(kEmpty) {}
    NodePtr(SingleLeaf* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kSingleLeaf) {}
    NodePtr(MultiLeaf* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kMultiLeaf) {}
    NodePtr(ForwardingNode* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kForwardingNode) {}
    NodePtr(BranchNode* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kBranchNode) {
      assert(ptr != nullptr);
    }

    Type getType() const { return Type(ptrAndType & 7u); }

    SingleLeaf* getSingleLeaf() const {
      assert(getType() == kSingleLeaf);
      return reinterpret_cast<SingleLeaf*>(ptrAndType & ~uintptr_t{7});
    }

    MultiLeaf* getMultiLeaf() const {
      assert(getType() == kMultiLeaf);
      return reinterpret_cast<MultiLeaf*>(ptrAndType & ~uintptr_t{7});
    }

    ForwardingNode* getForwardingNode() const {
      assert(getType() == kForwardingNode);
      return reinterpret_cast<ForwardingNode*>(ptrAndType & ~uintptr_t{7});
    }

    BranchNode* getBranchNode() const {
      assert(getType() == kBranchNode);
      return reinterpret_cast<BranchNode*>(ptrAndType & ~uintptr_t{7});
    }
  };

  // node type to skip ahead several levels to avoid needing to create a longer
  // chain of branching nodes between levels where the hash bits collide
  struct ForwardingNode {
    NodePtr targetNode;
    uint64_t hash;
    int numLevels;
  };

  struct BranchNode {
    uint64_t occupation;
    NodePtr child[1];
  };

  NodePtr root;

  bool insert_recurse(NodePtr* insertNode, uint64_t hash,
                      HighsHashTableEntry<K, V>& entry) {
    switch (insertNode->getType()) {
      case kEmpty: {
        *insertNode = new SingleLeaf(hash, std::move(entry));
        return true;
      }
      case kSingleLeaf: {
        SingleLeaf* leaf = insertNode->getSingleLeaf();
        if (leaf->hash == hash) {
          // check for existing key
          if (leaf->entry.key() == entry.key()) return false;
          // collision, turn single leaf into multi leaf
          MultiLeaf* newLeaf = new MultiLeaf(hash, std::move(entry));
          newLeaf->head.next = new ListNode(std::move(leaf->entry));
          *insertNode = newLeaf;
          // delete old leaf
          delete leaf;
          return true;
        }

        int numMatchingLevels = 0;
        uint64_t leafHash = leaf->hash;
        uint64_t currentHashBits = hash & 63;
        uint64_t leafHashBits = leaf->hash & 63;

        while (currentHashBits == leafHashBits) {
          ++numMatchingLevels;
          hash = hash >> 6;
          leafHash = leafHash >> 6;
          currentHashBits = hash & 63;
          leafHashBits = leafHash & 63;
        }

        BranchNode* newBranch =
            (BranchNode*)::operator new(sizeof(BranchNode) + sizeof(NodePtr));
        if (numMatchingLevels != 0) {
          ForwardingNode* forward = new ForwardingNode;
          *insertNode = forward;
          insertNode = &forward->targetNode;
          forward->numLevels = numMatchingLevels;
          forward->hash = leaf->hash;
        }

        leaf->hash = leafHash >> 6;
        int pos;
        newBranch->occupation =
            (uint64_t{1} << currentHashBits) | (uint64_t{1} << leafHashBits);
        pos = currentHashBits < leafHashBits;
        newBranch->child[pos] = nullptr;
        newBranch->child[1 - pos] = leaf;

        *insertNode = newBranch;
        insertNode = &newBranch->child[pos];
        hash = hash >> 6;
        break;
      }
      case kMultiLeaf: {
        MultiLeaf* leaf = insertNode->getMultiLeaf();
        if (leaf->hash == hash) {
          ListNode* iter = &leaf->head;

          while (true) {
            // check for existing key
            if (iter->entry.key() == entry.key()) return false;

            if (iter->next == nullptr) {
              // reached the end of the list and key is not duplicate, so insert
              iter->next = new ListNode(std::move(entry));
              return true;
            }
            iter = iter->next;
          }
        }

        int numMatchingLevels = 0;
        uint64_t leafHash = leaf->hash;
        uint64_t currentHashBits = hash & 63;
        uint64_t leafHashBits = leaf->hash & 63;

        while (currentHashBits == leafHashBits) {
          ++numMatchingLevels;
          hash = hash >> 6;
          leafHash = leafHash >> 6;
          currentHashBits = hash & 63;
          leafHashBits = leafHash & 63;
        }

        BranchNode* newBranch =
            (BranchNode*)::operator new(sizeof(BranchNode) + sizeof(NodePtr));
        if (numMatchingLevels != 0) {
          ForwardingNode* forward = new ForwardingNode;
          *insertNode = forward;
          insertNode = &forward->targetNode;
          forward->numLevels = numMatchingLevels;
          forward->hash = leaf->hash;
        }

        leaf->hash = leafHash >> 6;
        int pos;
        newBranch->occupation =
            (uint64_t{1} << currentHashBits) | (uint64_t{1} << leafHashBits);
        pos = currentHashBits < leafHashBits;
        newBranch->child[pos] = nullptr;
        newBranch->child[1 - pos] = leaf;

        *insertNode = newBranch;
        insertNode = &newBranch->child[pos];
        hash = hash >> 6;
        break;
      }
      case kForwardingNode: {
        ForwardingNode* forward = insertNode->getForwardingNode();

        int matchingLevels = 0;
        uint64_t unmatchedHash = forward->hash;
        do {
          if ((unmatchedHash & 63) != (hash & 63)) break;

          unmatchedHash = unmatchedHash >> 6;
          hash = hash >> 6;
          ++matchingLevels;
        } while (matchingLevels < forward->numLevels);

        if (matchingLevels == 0) {
          // we have an immediate mismatch, which means we create a new
          // branching node in front of the forwarding node and reduce the
          // forwarding nodes number of levels by one. If the forwarding node
          // only has a single level it's targetNode becomes a direct child of
          // the new branching node and the forwarding node is deleted.

          // the new branch node gets 2 children
          BranchNode* branch =
              (BranchNode*)::operator new(sizeof(BranchNode) + sizeof(NodePtr));
          int currentHashBits = hash & 63;
          int forwardHashBits = forward->hash & 63;
          branch->occupation = (uintptr_t{1} << currentHashBits) |
                               (uintptr_t{1} << forwardHashBits);
          int pos = currentHashBits < forwardHashBits;
          branch->child[pos] = nullptr;
          *insertNode = branch;
          insertNode = &branch->child[pos];
          hash = hash >> 6;
          --forward->numLevels;
          if (forward->numLevels == 0) {
            branch->child[1 - pos] = forward->targetNode;
            delete forward;
          } else {
            forward->hash = forward->hash >> 6;
            branch->child[1 - pos] = forward;
          }
        } else if (matchingLevels == forward->numLevels) {
          // we match all levels of the forwarding node, meaning we can simply
          // continue recursion and leave it untouched
          insertNode = &forward->targetNode;
        } else {
          // we need a new additional branching node and possibly a new
          // forwarding node to split the forwarded path with a branching in
          // between essentially replacing one of the previously forwarded
          // levels with a branching node
          int remainingLevels = forward->numLevels - matchingLevels - 1;
          forward->numLevels = matchingLevels;
          NodePtr oldTarget = forward->targetNode;

          BranchNode* branch =
              (BranchNode*)::operator new(sizeof(BranchNode) + sizeof(NodePtr));
          forward->targetNode = branch;

          int currentHashBits = hash & 63;
          int forwardHashBits = unmatchedHash & 63;
          branch->occupation = (uintptr_t{1} << currentHashBits) |
                               (uintptr_t{1} << forwardHashBits);

          int pos = currentHashBits < forwardHashBits;
          branch->child[pos] = nullptr;
          insertNode = &branch->child[pos];
          hash = hash >> 6;

          if (remainingLevels == 0) {
            branch->child[1 - pos] = oldTarget;
          } else {
            ForwardingNode* newForward = new ForwardingNode;
            branch->child[1 - pos] = newForward;
            newForward->numLevels = remainingLevels;
            newForward->targetNode = oldTarget;
            newForward->hash = unmatchedHash >> 6;
          }
        }

        break;
      }
      case kBranchNode: {
        BranchNode* branch = insertNode->getBranchNode();

        int bitPos = hash & 63;
        uint64_t bitMask = branch->occupation >> bitPos;
        int location = HighsHashHelpers::popcnt(bitMask);
        if (bitMask & 1) {
          --location;
        } else {
          int rightChilds =
              HighsHashHelpers::popcnt(branch->occupation << (63 - bitPos));
          assert(rightChilds + location ==
                 HighsHashHelpers::popcnt(branch->occupation));

          // sizeof(Branch) already contains the size for 1 pointer. So we just
          // need to add the left and right sizes up for the number of
          // additional pointers
          size_t leftSize = sizeof(BranchNode) + location * sizeof(NodePtr);
          size_t rightSize = rightChilds * sizeof(NodePtr);

          BranchNode* newBranch =
              (BranchNode*)::operator new(leftSize + rightSize);
          std::memcpy(newBranch, branch, leftSize - sizeof(NodePtr));
          std::memcpy(&newBranch->child[location + 1], &branch->child[location],
                      rightSize);
          newBranch->child[location] = nullptr;
          newBranch->occupation |= (uint64_t{1} << bitPos);
          *insertNode = newBranch;
          ::operator delete(branch);
          branch = newBranch;
        }

        insertNode = &branch->child[location];
        hash = hash >> 6;
      }
    }

    return insert_recurse(insertNode, hash, entry);
  }

  void erase_recurse(NodePtr* erase_node, uint64_t hash, const K& key) {
    switch (erase_node->getType()) {
      case kEmpty: {
        return;
      }
      case kSingleLeaf: {
        SingleLeaf* leaf = erase_node->getSingleLeaf();
        if (leaf->hash == hash) {
          // check for existing key
          if (leaf->entry.key() == key) {
            delete leaf;
            *erase_node = nullptr;
          }
        }
        return;
      }
      case kMultiLeaf: {
        MultiLeaf* leaf = erase_node->getMultiLeaf();
        if (leaf->hash == hash) {
          // check for existing key
          ListNode* iter = &leaf->head;

          do {
            ListNode* next = iter->next;
            if (iter->entry.key() == key) {
              if (next != nullptr) {
                *iter = std::move(*next);
                delete next;
                break;
              }
            }

            iter = next;
          } while (iter != nullptr);

          if (leaf->head.next == nullptr) {
            *erase_node =
                new SingleLeaf(leaf->hash, std::move(leaf->head.entry));
            delete leaf;
          }
        }

        break;
      }
      case kForwardingNode: {
        ForwardingNode* forward = erase_node->getForwardingNode();
        uint64_t mask = ((uint64_t{1} << (6 * forward->numLevels)) - 1);
        uint64_t forwardHash = forward->hash & mask;
        if (forwardHash != (hash & mask)) return;

        erase_recurse(&forward->targetNode, hash >> (6 * forward->numLevels),
                      key);
        switch (forward->targetNode.getType()) {
          case kEmpty:
            *erase_node = nullptr;
            delete forward;
            break;
          case kSingleLeaf: {
            *erase_node = forward->targetNode;
            SingleLeaf* leaf = erase_node->getSingleLeaf();
            leaf->hash = (leaf->hash << (6 * forward->numLevels)) | forwardHash;
            delete forward;
            break;
          }
          case kMultiLeaf: {
            *erase_node = forward->targetNode;
            MultiLeaf* leaf = erase_node->getMultiLeaf();
            leaf->hash = (leaf->hash << (6 * forward->numLevels)) | forwardHash;
            delete forward;
            break;
          }
          case kForwardingNode: {
            *erase_node = forward->targetNode;
            ForwardingNode* newForward = erase_node->getForwardingNode();
            newForward->hash =
                (newForward->hash << (6 * forward->numLevels)) | forwardHash;
            newForward->numLevels += forward->numLevels;
            delete forward;
          }
        }
        break;
      }
      case kBranchNode: {
        BranchNode* branch = erase_node->getBranchNode();

        int bitPos = hash & 63;
        uint64_t bitMask = branch->occupation >> bitPos;
        if ((bitMask & 1) == 0) return;

        int location = HighsHashHelpers::popcnt(bitMask) - 1;
        erase_recurse(&branch->child[location], hash >> 6, key);
        if (branch->child[location].getType() != kEmpty) return;

        branch->occupation = branch->occupation & (~(uint64_t{1} << bitPos));
        if (branch->occupation == 0) {
          *erase_node = nullptr;
          ::operator delete(branch);
          return;
        }

        int newSize = HighsHashHelpers::popcnt(branch->occupation);
        if (newSize == 1) {
          // there is one other entry so its index must be 1-location
          int pos = 1 - location;
          switch (branch->child[pos].getType()) {
            case kSingleLeaf: {
              SingleLeaf* leaf = branch->child[pos].getSingleLeaf();
              leaf->hash = (leaf->hash << 6) |
                           HighsHashHelpers::log2i(branch->occupation);
              *erase_node = leaf;
              break;
            }
            case kMultiLeaf: {
              MultiLeaf* leaf = branch->child[pos].getMultiLeaf();
              leaf->hash = (leaf->hash << 6) |
                           HighsHashHelpers::log2i(branch->occupation);
              *erase_node = leaf;
              break;
            }
            case kForwardingNode: {
              ForwardingNode* forward = branch->child[pos].getForwardingNode();
              ++forward->numLevels;
              forward->hash = (forward->hash << 6) |
                              HighsHashHelpers::log2i(branch->occupation);
              *erase_node = forward;
              break;
            }
            case kBranchNode: {
              ForwardingNode* forward = new ForwardingNode;
              forward->hash = HighsHashHelpers::log2i(branch->occupation);
              forward->targetNode = branch->child[pos];
              forward->numLevels = 1;
              *erase_node = forward;
              break;
            }
          }
        } else {
          BranchNode* compressedBranch = (BranchNode*)::operator new(
              sizeof(BranchNode) + (newSize - 1) * sizeof(NodePtr));
          *erase_node = compressedBranch;
          compressedBranch->occupation = branch->occupation;

          std::memcpy(&compressedBranch->child[0], &branch->child[0],
                      location * sizeof(NodePtr));
          std::memcpy(&compressedBranch->child[location],
                      &branch->child[location + 1],
                      (newSize - location) * sizeof(NodePtr));
        }

        ::operator delete(branch);
        break;
      }
    }
  }

  const ValueType* find_recurse(NodePtr node, uint64_t hash,
                                const K& key) const {
    switch (node.getType()) {
      case kEmpty:
        return nullptr;
      case kSingleLeaf:
        if (key == node.getSingleLeaf()->entry.key())
          return &node.getSingleLeaf()->entry.value();
        return nullptr;
      case kMultiLeaf: {
        MultiLeaf* leaf = node.getMultiLeaf();
        if (hash != leaf->hash) return nullptr;
        ListNode* iter = &leaf->head;
        do {
          if (iter->entry.key() == key) return &iter->entry.value();
          iter = iter->next;
        } while (iter != nullptr);
        return nullptr;
      }
      case kForwardingNode: {
        ForwardingNode* forward = node.getForwardingNode();

        uint64_t mask = (uint64_t{1} << (forward->numLevels * 6)) - 1;
        if (((hash ^ forward->hash) & mask) != 0) return nullptr;
        node = forward->targetNode;
        hash = hash >> (forward->numLevels * 6);
        break;
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int bitPos = hash & 63;
        uint64_t bitMask = branch->occupation >> bitPos;
        if ((bitMask & 1) == 0) return nullptr;
        int location = HighsHashHelpers::popcnt(bitMask) - 1;
        node = branch->child[location];
        hash = hash >> 6;
      }
    }

    return find_recurse(node, hash, key);
  }

  const HighsHashTableEntry<K,V>* find_common_recurse(NodePtr n1, NodePtr n2) const {
    if (n1.getType() > n2.getType()) std::swap(n1, n2);

    switch (n1.getType()) {
      case kEmpty:
        return nullptr;
      case kSingleLeaf: {
        SingleLeaf* leaf = n1.getSingleLeaf();
        if(find_recurse(n2, leaf->hash, leaf->entry.key()))
          return &leaf->entry;
        return nullptr;
      }
      case kMultiLeaf: {
        MultiLeaf* leaf = n1.getMultiLeaf();
        ListNode* iter = &leaf->head;
        do {
          if (find_recurse(n2, leaf->hash, iter->entry.key())) return &iter->entry;
          iter = iter->next;
        } while (iter != nullptr);
        return nullptr;
      }
      case kForwardingNode: {
        ForwardingNode* forward = n1.getForwardingNode();

        if (n2.getType() == kForwardingNode) {
          ForwardingNode* forward2 = n2.getForwardingNode();

          if (forward->numLevels > forward2->numLevels)
            std::swap(forward, forward2);

          uint64_t mask = (uint64_t{1} << (6 * forward->numLevels)) - 1;

          if ((forward->hash & mask) != (forward2->hash & mask)) return nullptr;

          if (forward2->numLevels == forward->numLevels) {
            return find_common_recurse(forward->targetNode,
                                       forward2->targetNode);
          }

          ForwardingNode tmpForward = *forward2;
          tmpForward.numLevels -= forward->numLevels;
          tmpForward.hash = tmpForward.hash >> (6 * forward->numLevels);

          return find_common_recurse(forward->targetNode, &tmpForward);
        }

        assert(n2.getType() == kBranchNode);

        BranchNode* branch = n2.getBranchNode();
        int bitPos = forward->hash & 63;
        uint64_t mask = branch->occupation >> bitPos;
        if ((mask & 1) == 0) return nullptr;
        int location = HighsHashHelpers::popcnt(mask) - 1;
        if (forward->numLevels == 1)
          return find_common_recurse(forward->targetNode,
                                     branch->child[location]);

        ForwardingNode tmpForward = *forward;
        --tmpForward.numLevels;
        tmpForward.hash = tmpForward.hash >> 6;

        return find_common_recurse(&tmpForward, branch->child[location]);
      }
      case kBranchNode: {
        BranchNode* branch1 = n1.getBranchNode();
        BranchNode* branch2 = n2.getBranchNode();

        uint64_t matchMask = branch1->occupation & branch2->occupation;

        while (matchMask) {
          int pos = HighsHashHelpers::log2i(matchMask);
          assert((branch1->occupation >> pos) & 1);
          assert((branch2->occupation >> pos) & 1);
          assert((matchMask >> pos) & 1);

          matchMask ^= (uint64_t{1} << pos);

          assert(((matchMask >> pos) & 1) == 0);

          int location1 =
              HighsHashHelpers::popcnt(branch1->occupation >> pos) - 1;
          int location2 =
              HighsHashHelpers::popcnt(branch2->occupation >> pos) - 1;

          const HighsHashTableEntry<K,V>* match = find_common_recurse(
              branch1->child[location1], branch2->child[location2]);
          if (match != nullptr) return match;
        }

        return nullptr;
      }
    }
  }

  void destroy_recurse(NodePtr node) {
    switch (node.getType()) {
      case kSingleLeaf:
        delete node.getSingleLeaf();
        break;
      case kMultiLeaf: {
        MultiLeaf* leaf = node.getMultiLeaf();
        delete leaf;
        ListNode* iter = leaf->head.next;
        while (iter != nullptr) {
          ListNode* next = iter->next;
          delete iter;
          iter = iter->next;
        }

        break;
      }
      case kForwardingNode: {
        ForwardingNode* forward = node.getForwardingNode();
        destroy_recurse(forward->targetNode);
        delete forward;
        break;
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = HighsHashHelpers::popcnt(branch->occupation);

        for (int i = 0; i < size; ++i) destroy_recurse(branch->child[i]);

        ::operator delete(branch);
      }
    }
  }

  NodePtr copy_recurse(NodePtr node) {
    switch (node.getType()) {
      case kSingleLeaf:
        return new SingleLeaf(*node.getSingleLeaf());
      case kMultiLeaf: {
        MultiLeaf* leaf = node.getMultiLeaf();

        MultiLeaf* copyLeaf = new MultiLeaf(*leaf);

        ListNode* iter = &leaf->head;
        ListNode* copyIter = &copyLeaf->head;
        do {
          copyIter->next = new ListNode(*iter->next);
          iter = iter->next;
          copyIter = copyIter->next;
        } while (iter->next != nullptr);

        return copyLeaf;
      }
      case kForwardingNode: {
        ForwardingNode* forward = new ForwardingNode(*node.getForwardingNode());
        forward->targetNode = copy_recurse(forward->targetNode);
        return forward;
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = HighsHashHelpers::popcnt(branch->occupation);
        BranchNode* newBranch = (BranchNode*)::operator new(
            sizeof(BranchNode) + (size - 1) * sizeof(NodePtr));
        newBranch->occupation = branch->occupation;
        for (int i = 0; i < size; ++i)
          newBranch->child[i] = copy_recurse(branch->child[i]);

        return newBranch;
      }
    }
  }

  template <typename F>
  bool for_each_recurse(NodePtr node, F&& f) const {
    switch (node.getType()) {
      case kEmpty:
        //abort();
        break;
      case kSingleLeaf: {
        SingleLeaf* leaf = node.getSingleLeaf();
        assert(leaf != nullptr);
        return f(leaf->entry);
      }
      case kMultiLeaf: {
        MultiLeaf* leaf = node.getMultiLeaf();
        ListNode* iter = &leaf->head;
        do {
          if (f(iter->entry)) return true;
          iter = iter->next;
        } while (iter != nullptr);
        break;
      }
      case kForwardingNode: {
        ForwardingNode* forward = node.getForwardingNode();
        return for_each_recurse(forward->targetNode, f);
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = HighsHashHelpers::popcnt(branch->occupation);

        for (int i = 0; i < size; ++i)
          if (for_each_recurse(branch->child[i], f)) return true;
      }
    }

    return false;
  }

 public:
  template <typename... Args>
  bool insert(Args&&... args) {
    HighsHashTableEntry<K, V> entry(std::forward<Args>(args)...);
    uint64_t hash = HighsHashHelpers::hash(entry.key());
    return insert_recurse(&root, hash, entry);
  }

  void erase(const K& key) {
    uint64_t hash = HighsHashHelpers::hash(key);

    erase_recurse(&root, hash, key);
  }

  bool contains(const K& key) const {
    uint64_t hash = HighsHashHelpers::hash(key);

    return find_recurse(root, hash, key) != nullptr;
  }

  const ValueType* find(const K& key) const {
    uint64_t hash = HighsHashHelpers::hash(key);

    return find_recurse(root, hash, key);
  }

  const HighsHashTableEntry<K,V>* find_common(const HighsHashTree<K, V>& other) const {
    return find_common_recurse(root, other.root);
  }

  bool empty() const { return root.getType() == kEmpty; }

  void clear() {
    destroy_recurse(root);
    root = nullptr;
  }

  template <typename F>
  bool for_each(F&& f) const {
    return for_each_recurse(root, f);
  }

  HighsHashTree() = default;

  HighsHashTree(HighsHashTree&& other) : root(other.root) {
    other.root = nullptr;
  }

  HighsHashTree(const HighsHashTree& other) : root(copy_recurse(other.root)) {}

  HighsHashTree& operator=(HighsHashTree&& other) {
    destroy_recurse(root);
    root = other.root;
    other.root = nullptr;
    return *this;
  }

  HighsHashTree& operator=(const HighsHashTree& other) {
    destroy_recurse(root);
    root = copy_recurse(other.root);
    return *this;
  }

  ~HighsHashTree() { destroy_recurse(root); }
};

#endif
