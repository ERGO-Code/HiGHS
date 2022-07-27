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

using std::memcpy;
using std::memmove;

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

  enum Constants {
    kBitsPerLevel = 6,
    kBranchFactor = 1 << kBitsPerLevel,
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
    int targetLevel;
  };

  struct BranchNode {
    uint64_t occupation;
    NodePtr child[1];

    void set_occupation(uint8_t pos) { occupation |= uint64_t{1} << (pos); }

    void flip_occupation(uint8_t pos) { occupation ^= uint64_t{1} << (pos); }

    bool test_occupation(uint8_t pos) const {
      return occupation & (uint64_t{1} << pos);
    }

    int get_child_position(uint8_t pos) const {
      return HighsHashHelpers::popcnt(occupation >> pos);
    }

    int get_num_children_after(uint8_t pos) const {
      return HighsHashHelpers::popcnt(occupation << (63 - (pos)));
    }

    int get_num_children() const {
      return HighsHashHelpers::popcnt(occupation);
    }
  };

  // allocate branch nodes in multiples of 64 bytes to reduce allocator stress
  // with different sizes and reduce reallocations of nodes
  static constexpr size_t getBranchNodeSize(int numChilds) {
    return (sizeof(BranchNode) + size_t(numChilds - 1) * sizeof(NodePtr) + 63) &
           ~63;
  };

  static ForwardingNode* createForwardingNode(uint64_t hash, int targetLevel) {
    ForwardingNode* forward =
        (ForwardingNode*)::operator new(getBranchNodeSize(2));
    forward->hash = hash;
    forward->targetLevel = targetLevel;
    return forward;
  }
  static BranchNode* createBranchingNode(int numChilds) {
    BranchNode* branch =
        (BranchNode*)::operator new(getBranchNodeSize(numChilds));
    branch->occupation = 0;
    return branch;
  }

  static BranchNode* convertForwardingNodeToBinaryBranch(void* forward) {
    BranchNode* branch = (BranchNode*)forward;
    branch->occupation = 0;
    return branch;
  }

  static ForwardingNode* convertUnaryBranchToForwardingNode(
      void* branch, NodePtr singleChild, uint64_t hash, int branchNodeLevel) {
    ForwardingNode* forward = (ForwardingNode*)branch;

    forward->hash = hash;
    forward->targetNode = singleChild;
    forward->targetLevel = branchNodeLevel;

    return forward;
  }

  static void destroyInnerNode(void* innerNode) {
    ::operator delete(innerNode);
  }

  static BranchNode* addChildToBranchNode(BranchNode* branch, uint8_t hashValue,
                                          int location) {
    int rightChilds = branch->get_num_children_after(hashValue);
    assert(rightChilds + location == branch->get_num_children());

    size_t newSize = getBranchNodeSize(location + rightChilds + 1);
    size_t rightSize = rightChilds * sizeof(NodePtr);

    if (newSize == getBranchNodeSize(location + rightChilds)) {
      memmove(&branch->child[location + 1], &branch->child[location],
              rightSize);

      return branch;
    }

    BranchNode* newBranch = (BranchNode*)::operator new(newSize);
    // sizeof(Branch) already contains the size for 1 pointer. So we just
    // need to add the left and right sizes up for the number of
    // additional pointers
    size_t leftSize = sizeof(BranchNode) + (location - 1) * sizeof(NodePtr);

    memcpy(newBranch, branch, leftSize);
    memcpy(&newBranch->child[location + 1], &branch->child[location],
           rightSize);

    destroyInnerNode(branch);

    return newBranch;
  }

  static NodePtr removeChildFromBranchNode(BranchNode* branch, int location,
                                           uint64_t hash, int hashPos) {
    NodePtr newNode;
    int newNumChild = branch->get_num_children();

    if (newNumChild == 1) {
      // there is one other entry so its index must be 1-location
      int pos = 1 - location;
      switch (branch->child[pos].getType()) {
        case kSingleLeaf:
        case kMultiLeaf:
        case kForwardingNode:
          newNode = branch->child[pos];
          destroyInnerNode(branch);
          break;
        case kBranchNode: {
          uint8_t childHash = HighsHashHelpers::log2i(branch->occupation);

          assert(branch->get_child_position(childHash) == 1);
          assert(branch->test_occupation(childHash));

          ForwardingNode* forward = convertUnaryBranchToForwardingNode(
              branch, branch->child[pos], hash, hashPos + 1);

          set_hash_chunk(forward->hash, childHash, hashPos);
          newNode = forward;
        }
      }
    } else {
      size_t newSize = getBranchNodeSize(newNumChild);
      size_t rightSize = (newNumChild - location) * sizeof(NodePtr);
      if (newSize == getBranchNodeSize(newNumChild + 1)) {
        // allocated size class is the same, so we do not allocate a new node
        memmove(&branch->child[location], &branch->child[location + 1],
                rightSize);
        newNode = branch;
      } else {
        // allocated size class changed, so we allocate a smaller branch node
        BranchNode* compressedBranch = (BranchNode*)::operator new(newSize);
        newNode = compressedBranch;

        size_t leftSize =
            offsetof(BranchNode, child) + location * sizeof(NodePtr);
        memcpy(compressedBranch, branch, leftSize);
        memcpy(&compressedBranch->child[location], &branch->child[location + 1],
               rightSize);

        destroyInnerNode(branch);
      }
    }

    return newNode;
  }

  NodePtr root;

  static uint8_t get_hash_chunk(uint64_t hash, int pos) {
    return (hash >> (pos * kBitsPerLevel)) & (kBranchFactor - 1);
  }

  static void set_hash_chunk(uint64_t& hash, uint64_t chunk, int chunkPos) {
    const int shiftAmount = chunkPos * kBitsPerLevel;
    chunk ^= (hash >> shiftAmount) & (kBranchFactor - 1);
    hash ^= chunk << shiftAmount;
  }

  static uint64_t get_first_n_hash_chunks(uint64_t hash, int n) {
    return hash & ((uint64_t{1} << (n * kBitsPerLevel)) - 1);
  }

  bool insert_recurse(NodePtr* insertNode, uint64_t hash, int hashPos,
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

        int diffPos = hashPos;
        while (get_hash_chunk(hash, diffPos) ==
               get_hash_chunk(leaf->hash, diffPos))
          ++diffPos;

        if (diffPos > hashPos) {
          ForwardingNode* forward = createForwardingNode(leaf->hash, diffPos);
          *insertNode = forward;
          insertNode = &forward->targetNode;
          hashPos = diffPos;
        }

        BranchNode* newBranch = createBranchingNode(2);
        newBranch->set_occupation(get_hash_chunk(hash, hashPos));
        newBranch->set_occupation(get_hash_chunk(leaf->hash, hashPos));

        int pos =
            get_hash_chunk(hash, hashPos) < get_hash_chunk(leaf->hash, hashPos);
        newBranch->child[pos] = nullptr;
        newBranch->child[1 - pos] = leaf;

        *insertNode = newBranch;
        insertNode = &newBranch->child[pos];
        ++hashPos;
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

        int diffPos = hashPos;
        while (get_hash_chunk(hash, diffPos) ==
               get_hash_chunk(leaf->hash, diffPos))
          ++diffPos;

        if (diffPos > hashPos) {
          ForwardingNode* forward = createForwardingNode(leaf->hash, diffPos);
          *insertNode = forward;
          insertNode = &forward->targetNode;
          hashPos = diffPos;
        }

        BranchNode* newBranch = createBranchingNode(2);
        newBranch->set_occupation(get_hash_chunk(hash, hashPos));
        newBranch->set_occupation(get_hash_chunk(leaf->hash, hashPos));

        int pos =
            get_hash_chunk(hash, hashPos) < get_hash_chunk(leaf->hash, hashPos);
        newBranch->child[pos] = nullptr;
        newBranch->child[1 - pos] = leaf;

        *insertNode = newBranch;
        insertNode = &newBranch->child[pos];
        ++hashPos;
        break;
      }
      case kForwardingNode: {
        ForwardingNode* forward = insertNode->getForwardingNode();

        int diffPos = hashPos;
        do {
          if (get_hash_chunk(hash, diffPos) !=
              get_hash_chunk(forward->hash, diffPos))
            break;
          ++diffPos;
        } while (diffPos < forward->targetLevel);

        if (diffPos == hashPos) {
          // we have an immediate mismatch, which means we create a new
          // branching node in front of the forwarding node and reduce the
          // forwarding nodes number of levels by one. If the forwarding node
          // only has a single level it's targetNode becomes a direct child of
          // the new branching node and the forwarding node is deleted.

          BranchNode* branch;
          uint8_t forwardDiffHash = get_hash_chunk(forward->hash, diffPos);
          int pos = get_hash_chunk(hash, diffPos) < forwardDiffHash;

          if (forward->targetLevel - diffPos == 1) {
            NodePtr targetNode = forward->targetNode;
            // forwarding node can be converted to branch node without new
            // allocation
            branch = convertForwardingNodeToBinaryBranch(forward);
            branch->child[1 - pos] = targetNode;
          } else {
            branch = createBranchingNode(2);
            branch->child[1 - pos] = forward;
          }

          branch->set_occupation(forwardDiffHash);
          branch->set_occupation(get_hash_chunk(hash, diffPos));
          branch->child[pos] = nullptr;
          *insertNode = branch;
          insertNode = &branch->child[pos];
          hashPos = diffPos + 1;
        } else if (diffPos == forward->targetLevel) {
          // we match all levels of the forwarding node, meaning we can simply
          // continue recursion and leave it untouched
          hashPos = diffPos;
          insertNode = &forward->targetNode;
        } else {
          // we need a new additional branching node and possibly a new
          // forwarding node to split the forwarded path with a branching in
          // between essentially replacing one of the previously forwarded
          // levels with a branching node
          int oldSplit = forward->targetLevel;
          forward->targetLevel = diffPos;
          NodePtr oldTarget = forward->targetNode;

          BranchNode* branch = createBranchingNode(2);
          forward->targetNode = branch;
          branch->set_occupation(get_hash_chunk(hash, diffPos));
          branch->set_occupation(get_hash_chunk(forward->hash, diffPos));

          int pos = get_hash_chunk(hash, diffPos) <
                    get_hash_chunk(forward->hash, diffPos);
          branch->child[pos] = nullptr;
          insertNode = &branch->child[pos];
          hashPos = diffPos + 1;

          if (hashPos == oldSplit) {
            branch->child[1 - pos] = oldTarget;
          } else {
            ForwardingNode* newForward =
                createForwardingNode(forward->hash, oldSplit);
            branch->child[1 - pos] = newForward;
            newForward->targetNode = oldTarget;
          }
        }

        break;
      }
      case kBranchNode: {
        BranchNode* branch = insertNode->getBranchNode();

        int location =
            branch->get_child_position(get_hash_chunk(hash, hashPos));

        if (branch->test_occupation(get_hash_chunk(hash, hashPos))) {
          --location;
        } else {
          branch = addChildToBranchNode(branch, get_hash_chunk(hash, hashPos),
                                        location);

          branch->child[location] = nullptr;
          branch->set_occupation(get_hash_chunk(hash, hashPos));
        }

        *insertNode = branch;
        insertNode = &branch->child[location];
        ++hashPos;
      }
    }

    return insert_recurse(insertNode, hash, hashPos, entry);
  }

  void erase_recurse(NodePtr* erase_node, uint64_t hash, int hashPos,
                     const K& key) {
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
        if (get_first_n_hash_chunks(forward->hash, forward->targetLevel) !=
            get_first_n_hash_chunks(hash, forward->targetLevel))
          return;

        erase_recurse(&forward->targetNode, hash, forward->targetLevel, key);
        switch (forward->targetNode.getType()) {
          case kEmpty:
            *erase_node = nullptr;
            destroyInnerNode(forward);
            break;
          case kSingleLeaf:
          case kMultiLeaf:
          case kForwardingNode:
            *erase_node = forward->targetNode;
            destroyInnerNode(forward);
        }
        break;
      }
      case kBranchNode: {
        BranchNode* branch = erase_node->getBranchNode();

        if (!branch->test_occupation(get_hash_chunk(hash, hashPos))) return;

        int location =
            branch->get_child_position(get_hash_chunk(hash, hashPos)) - 1;
        erase_recurse(&branch->child[location], hash, hashPos + 1, key);

        if (branch->child[location].getType() != kEmpty) return;

        branch->flip_occupation(get_hash_chunk(hash, hashPos));

        *erase_node =
            removeChildFromBranchNode(branch, location, hash, hashPos);
        break;
      }
    }
  }

  const ValueType* find_recurse(NodePtr node, uint64_t hash, int hashPos,
                                const K& key) const {
    int startPos = hashPos;
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
        if (get_first_n_hash_chunks(hash, forward->targetLevel) !=
            get_first_n_hash_chunks(forward->hash, forward->targetLevel))
          return nullptr;

        node = forward->targetNode;
        assert(forward->targetLevel > hashPos);
        hashPos = forward->targetLevel;
        break;
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        if (!branch->test_occupation(get_hash_chunk(hash, hashPos)))
          return nullptr;
        int location =
            branch->get_child_position(get_hash_chunk(hash, hashPos)) - 1;
        node = branch->child[location];
        ++hashPos;
      }
    }

    assert(hashPos > startPos);

    return find_recurse(node, hash, hashPos, key);
  }

  const HighsHashTableEntry<K, V>* find_common_recurse(NodePtr n1, NodePtr n2,
                                                       int hashPos) const {
    if (n1.getType() > n2.getType()) std::swap(n1, n2);

    switch (n1.getType()) {
      case kEmpty:
        return nullptr;
      case kSingleLeaf: {
        SingleLeaf* leaf = n1.getSingleLeaf();
        if (find_recurse(n2, leaf->hash, hashPos, leaf->entry.key()))
          return &leaf->entry;
        return nullptr;
      }
      case kMultiLeaf: {
        MultiLeaf* leaf = n1.getMultiLeaf();
        ListNode* iter = &leaf->head;
        do {
          if (find_recurse(n2, leaf->hash, hashPos, iter->entry.key()))
            return &iter->entry;
          iter = iter->next;
        } while (iter != nullptr);
        return nullptr;
      }
      case kForwardingNode: {
        ForwardingNode* forward = n1.getForwardingNode();

        if (n2.getType() == kForwardingNode) {
          ForwardingNode* forward2 = n2.getForwardingNode();

          if (forward->targetLevel > forward2->targetLevel)
            std::swap(forward, forward2);

          if (get_first_n_hash_chunks(forward->hash, forward->targetLevel) !=
              get_first_n_hash_chunks(forward2->hash, forward->targetLevel))
            return nullptr;

          return find_common_recurse(
              forward->targetNode,
              forward2->targetLevel == forward->targetLevel
                  ? forward2->targetNode
                  : n2,
              forward->targetLevel);
        }

        assert(n2.getType() == kBranchNode);

        BranchNode* branch = n2.getBranchNode();
        if (!branch->test_occupation(get_hash_chunk(forward->hash, hashPos)))
          return nullptr;
        int location =
            branch->get_child_position(get_hash_chunk(forward->hash, hashPos)) -
            1;

        ++hashPos;
        return find_common_recurse(
            forward->targetLevel == hashPos ? forward->targetNode : n1,
            branch->child[location], hashPos);
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

          int location1 = branch1->get_child_position(pos) - 1;
          int location2 = branch2->get_child_position(pos) - 1;

          const HighsHashTableEntry<K, V>* match =
              find_common_recurse(branch1->child[location1],
                                  branch2->child[location2], hashPos + 1);
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
        destroyInnerNode(forward);
        break;
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = branch->get_num_children();

        for (int i = 0; i < size; ++i) destroy_recurse(branch->child[i]);

        destroyInnerNode(branch);
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
        int size = branch->get_num_children();
        BranchNode* newBranch =
            (BranchNode*)::operator new(getBranchNodeSize(size));
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
        int size = branch->get_num_children();

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
    return insert_recurse(&root, hash, 0, entry);
  }

  void erase(const K& key) {
    uint64_t hash = HighsHashHelpers::hash(key);

    erase_recurse(&root, hash, 0, key);
  }

  bool contains(const K& key) const {
    uint64_t hash = HighsHashHelpers::hash(key);
    return find_recurse(root, hash, 0, key) != nullptr;
  }

  const ValueType* find(const K& key) const {
    uint64_t hash = HighsHashHelpers::hash(key);

    return find_recurse(root, hash, 0, key);
  }

  const HighsHashTableEntry<K, V>* find_common(
      const HighsHashTree<K, V>& other) const {
    return find_common_recurse(root, other.root, 0);
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
