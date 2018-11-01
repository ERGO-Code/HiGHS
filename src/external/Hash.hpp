/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef EXTERNAL_HASH_HPP_
#define EXTERNAL_HASH_HPP_

#include "tessil/hopscotch_map.h"
#include "tessil/hopscotch_set.h"
#include <cstdint>

#if 1 // TODO check performance of applying the multiplicative on top of
      // std::hash, std::hash often uses just the indentity

template <typename T, typename H>
struct Hasher
{
   size_t
   operator()( const T& x ) const
   {
      uint64_t y = H{}( x );
      return size_t( ( y * uint64_t( 0x9e3779b97f4a7c15 ) ) >> 32 );
   }
};

template <typename K, typename V, typename H = std::hash<K>>
using HashMap = tsl::hopscotch_map<K, V, Hasher<K, H>>;

template <typename T, typename H = std::hash<T>>
using HashSet = tsl::hopscotch_set<T, Hasher<T, H>>;

#else

template <typename K, typename V, typename H = std::hash<K>>
using HashMap = tsl::hopscotch_map<K, V, H>;

template <typename T, typename H = std::hash<T>>
using HashSet = tsl::hopscotch_set<T, H>;

#endif

#endif