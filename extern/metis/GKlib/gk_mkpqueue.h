/*!
\file  gk_mkpqueue.h
\brief Templates for priority queues

\date   Started 4/09/07
\author George
\version\verbatim $Id: gk_mkpqueue.h 21742 2018-01-26 16:59:15Z karypis $ \endverbatim
*/


#ifndef _GK_MKPQUEUE_H
#define _GK_MKPQUEUE_H


#define GK_MKPQUEUE(FPRFX, PQT, KVT, KT, VT, KVMALLOC, KMAX, KEY_LT)\
/*************************************************************************/\
/*! This function creates and initializes a priority queue */\
/**************************************************************************/\
PQT *FPRFX ## Create(size_t maxnodes)\
{\
  PQT *queue; \
\
  queue = (PQT *)malloc(sizeof(PQT));\
  FPRFX ## Init(queue, maxnodes);\
\
  return queue;\
}\
\
\
/*************************************************************************/\
/*! This function initializes the data structures of the priority queue */\
/**************************************************************************/\
void FPRFX ## Init(PQT *queue, size_t maxnodes)\
{\
  queue->nnodes = 0;\
  queue->maxnodes = maxnodes;\
\
  queue->heap    = KVMALLOC(maxnodes);\
  queue->locator = (gk_idx_t*)malloc(maxnodes*sizeof(gk_idx_t));\
  for (size_t i = 0; i < maxnodes; ++i)\
    queue->locator[i] = -1;\
}\
\
\
/*************************************************************************/\
/*! This function resets the priority queue */\
/**************************************************************************/\
void FPRFX ## Reset(PQT *queue)\
{\
  ssize_t i;\
  ssize_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  for (i=queue->nnodes-1; i>=0; i--)\
    locator[heap[i].val] = -1;\
  queue->nnodes = 0;\
}\
\
\
/*************************************************************************/\
/*! This function frees the internal datastructures of the priority queue */\
/**************************************************************************/\
void FPRFX ## Free(PQT *queue)\
{\
  if (queue == NULL) return;\
  gk_free((void **)&queue->heap);\
  gk_free((void **)&queue->locator);\
  queue->maxnodes = 0;\
}\
\
\
/*************************************************************************/\
/*! This function frees the internal datastructures of the priority queue \
    and the queue itself */\
/**************************************************************************/\
void FPRFX ## Destroy(PQT *queue)\
{\
  if (queue == NULL) return;\
  FPRFX ## Free(queue);\
  gk_free((void **)&queue);\
}\
\
\
/*************************************************************************/\
/*! This function returns the length of the queue */\
/**************************************************************************/\
size_t FPRFX ## Length(PQT *queue)\
{\
  return queue->nnodes;\
}\
\
\
/*************************************************************************/\
/*! This function adds an item in the priority queue */\
/**************************************************************************/\
int FPRFX ## Insert(PQT *queue, VT node, KT key)\
{\
  ssize_t i, j;\
  ssize_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  i = queue->nnodes++;\
  while (i > 0) {\
    j = (i-1)>>1;\
    if (KEY_LT(key, heap[j].key)) {\
      heap[i] = heap[j];\
      locator[heap[i].val] = i;\
      i = j;\
    }\
    else\
      break;\
  }\
  heap[i].key   = key;\
  heap[i].val   = node;\
  locator[node] = i;\
\
  return 0;\
}\
\
\
/*************************************************************************/\
/*! This function deletes an item from the priority queue */\
/**************************************************************************/\
int FPRFX ## Delete(PQT *queue, VT node)\
{\
  ssize_t i, j;\
  size_t nnodes;\
  KT newkey, oldkey;\
  ssize_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  i = locator[node];\
  locator[node] = -1;\
\
  if (--queue->nnodes > 0 && heap[queue->nnodes].val != node) {\
    node   = heap[queue->nnodes].val;\
    newkey = heap[queue->nnodes].key;\
    oldkey = heap[i].key;\
\
    if (KEY_LT(newkey, oldkey)) { /* Filter-up */\
      while (i > 0) {\
        j = (i-1)>>1;\
        if (KEY_LT(newkey, heap[j].key)) {\
          heap[i] = heap[j];\
          locator[heap[i].val] = i;\
          i = j;\
        }\
        else\
          break;\
      }\
    }\
    else { /* Filter down */\
      nnodes = queue->nnodes;\
      while ((j=(i<<1)+1) < nnodes) {\
        if (KEY_LT(heap[j].key, newkey)) {\
          if (j+1 < nnodes && KEY_LT(heap[j+1].key, heap[j].key))\
            j++;\
          heap[i] = heap[j];\
          locator[heap[i].val] = i;\
          i = j;\
        }\
        else if (j+1 < nnodes && KEY_LT(heap[j+1].key, newkey)) {\
          j++;\
          heap[i] = heap[j];\
          locator[heap[i].val] = i;\
          i = j;\
        }\
        else\
          break;\
      }\
    }\
\
    heap[i].key   = newkey;\
    heap[i].val   = node;\
    locator[node] = i;\
  }\
\
  return 0;\
}\
\
\
/*************************************************************************/\
/*! This function updates the key values associated for a particular item */ \
/**************************************************************************/\
void FPRFX ## Update(PQT *queue, VT node, KT newkey)\
{\
  ssize_t i, j;\
  size_t nnodes;\
  KT oldkey;\
  ssize_t *locator=queue->locator;\
  KVT *heap=queue->heap;\
\
  oldkey = heap[locator[node]].key;\
  if (!KEY_LT(newkey, oldkey) && !KEY_LT(oldkey, newkey)) return;\
\
  i = locator[node];\
\
  if (KEY_LT(newkey, oldkey)) { /* Filter-up */\
    while (i > 0) {\
      j = (i-1)>>1;\
      if (KEY_LT(newkey, heap[j].key)) {\
        heap[i] = heap[j];\
        locator[heap[i].val] = i;\
        i = j;\
      }\
      else\
        break;\
    }\
  }\
  else { /* Filter down */\
    nnodes = queue->nnodes;\
    while ((j=(i<<1)+1) < nnodes) {\
      if (KEY_LT(heap[j].key, newkey)) {\
        if (j+1 < nnodes && KEY_LT(heap[j+1].key, heap[j].key))\
          j++;\
        heap[i] = heap[j];\
        locator[heap[i].val] = i;\
        i = j;\
      }\
      else if (j+1 < nnodes && KEY_LT(heap[j+1].key, newkey)) {\
        j++;\
        heap[i] = heap[j];\
        locator[heap[i].val] = i;\
        i = j;\
      }\
      else\
        break;\
    }\
  }\
\
  heap[i].key   = newkey;\
  heap[i].val   = node;\
  locator[node] = i;\
\
  return;\
}\
\
\
/*************************************************************************/\
/*! This function returns the item at the top of the queue and removes\
    it from the priority queue */\
/**************************************************************************/\
VT FPRFX ## GetTop(PQT *queue)\
{\
  ssize_t i, j;\
  ssize_t *locator;\
  KVT *heap;\
  VT vtx, node;\
  KT key;\
\
  if (queue->nnodes == 0)\
    return -1;\
\
  queue->nnodes--;\
\
  heap    = queue->heap;\
  locator = queue->locator;\
\
  vtx = heap[0].val;\
  locator[vtx] = -1;\
\
  if ((i = queue->nnodes) > 0) {\
    key  = heap[i].key;\
    node = heap[i].val;\
    i = 0;\
    while ((j=2*i+1) < queue->nnodes) {\
      if (KEY_LT(heap[j].key, key)) {\
        if (j+1 < queue->nnodes && KEY_LT(heap[j+1].key, heap[j].key))\
          j = j+1;\
        heap[i] = heap[j];\
        locator[heap[i].val] = i;\
        i = j;\
      }\
      else if (j+1 < queue->nnodes && KEY_LT(heap[j+1].key, key)) {\
        j = j+1;\
        heap[i] = heap[j];\
        locator[heap[i].val] = i;\
        i = j;\
      }\
      else\
        break;\
    }\
\
    heap[i].key   = key;\
    heap[i].val   = node;\
    locator[node] = i;\
  }\
\
  return vtx;\
}\
\
\
/*************************************************************************/\
/*! This function returns the item at the top of the queue. The item is not\
    deleted from the queue. */\
/**************************************************************************/\
VT FPRFX ## SeeTopVal(PQT *queue)\
{\
  return (queue->nnodes == 0 ? -1 : queue->heap[0].val);\
}\
\
\
/*************************************************************************/\
/*! This function returns the key of the top item. The item is not\
    deleted from the queue. */\
/**************************************************************************/\
KT FPRFX ## SeeTopKey(PQT *queue)\
{\
  return (queue->nnodes == 0 ? KMAX : queue->heap[0].key);\
}\


#define GK_MKPQUEUE_PROTO(FPRFX, PQT, KT, VT)\
  PQT *  FPRFX ## Create(size_t maxnodes);\
  void   FPRFX ## Init(PQT *queue, size_t maxnodes);\
  void   FPRFX ## Reset(PQT *queue);\
  void   FPRFX ## Free(PQT *queue);\
  void   FPRFX ## Destroy(PQT *queue);\
  size_t FPRFX ## Length(PQT *queue);\
  int    FPRFX ## Insert(PQT *queue, VT node, KT key);\
  int    FPRFX ## Delete(PQT *queue, VT node);\
  void   FPRFX ## Update(PQT *queue, VT node, KT newkey);\
  VT     FPRFX ## GetTop(PQT *queue);\
  VT     FPRFX ## SeeTopVal(PQT *queue);\
  KT     FPRFX ## SeeTopKey(PQT *queue);\



#endif
