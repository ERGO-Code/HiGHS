#ifndef HIGHS_METIS_MANGLE_H
#define HIGHS_METIS_MANGLE_H

/* balance.c */
#define Balance2Way                         Highs_metis_Balance2Way
#define Bnd2WayBalance                      Highs_metis_Bnd2WayBalance
#define General2WayBalance                  Highs_metis_General2WayBalance
#define McGeneral2WayBalance                Highs_metis_McGeneral2WayBalance

/* bucketsort.c */
#define BucketSortKeysInc                   Highs_metis_BucketSortKeysInc

/* coarsen.c */
#define CoarsenGraph                        Highs_metis_CoarsenGraph
#define CoarsenGraphNlevels                 Highs_metis_CoarsenGraphNlevels
#define Match_RM                            Highs_metis_Match_RM
#define Match_SHEM                          Highs_metis_Match_SHEM
#define Match_2Hop                          Highs_metis_Match_2Hop
#define Match_2HopAny                       Highs_metis_Match_2HopAny
#define Match_2HopAll                       Highs_metis_Match_2HopAll
#define PrintCGraphStats                    Highs_metis_PrintCGraphStats
#define CreateCoarseGraph                   Highs_metis_CreateCoarseGraph
#define SetupCoarseGraph                    Highs_metis_SetupCoarseGraph
#define ReAdjustMemory                      Highs_metis_ReAdjustMemory

/* compress.c */
#define CompressGraph                       Highs_metis_CompressGraph
#define PruneGraph                          Highs_metis_PruneGraph

/* contig.c */
#define FindSepInducedComponents            Highs_metis_FindSepInducedComponents

/* fm.c */
#define FM_2WayRefine                       Highs_metis_FM_2WayRefine
#define FM_2WayCutRefine                    Highs_metis_FM_2WayCutRefine
#define FM_Mc2WayCutRefine                  Highs_metis_FM_Mc2WayCutRefine
#define SelectQueue                         Highs_metis_SelectQueue
#define Print2WayRefineStats                Highs_metis_Print2WayRefineStats

/* graph.c */
#define SetupGraph                          Highs_metis_SetupGraph
#define SetupGraph_tvwgt                    Highs_metis_SetupGraph_tvwgt
#define SetupGraph_label                    Highs_metis_SetupGraph_label
#define SetupSplitGraph                     Highs_metis_SetupSplitGraph
#define CreateGraph                         Highs_metis_CreateGraph
#define InitGraph                           Highs_metis_InitGraph
#define FreeSData                           Highs_metis_FreeSData
#define FreeRData                           Highs_metis_FreeRData
#define FreeGraph                           Highs_metis_FreeGraph

/* initpart.c */
#define InitSeparator                       Highs_metis_InitSeparator
#define RandomBisection                     Highs_metis_RandomBisection
#define GrowBisection                       Highs_metis_GrowBisection
#define GrowBisectionNode                   Highs_metis_GrowBisectionNode

/* mcutil.c */
#define ivecle                              Highs_metis_ivecle
#define ivecaxpylez                         Highs_metis_ivecaxpylez
#define BetterVBalance                      Highs_metis_BetterVBalance
#define BetterBalance2Way                   Highs_metis_BetterBalance2Way
#define ComputeLoadImbalance                Highs_metis_ComputeLoadImbalance
#define ComputeLoadImbalanceDiff            Highs_metis_ComputeLoadImbalanceDiff
#define ComputeLoadImbalanceDiffVec         Highs_metis_ComputeLoadImbalanceDiffVec

/* mmd.c */
#define genmmd                              Highs_metis_genmmd
#define mmdelm                              Highs_metis_mmdelm
#define mmdint                              Highs_metis_mmdint
#define mmdnum                              Highs_metis_mmdnum
#define mmdupd                              Highs_metis_mmdupd

/* ometis.c */
#define MlevelNestedDissection              Highs_metis_MlevelNestedDissection
#define MlevelNestedDissectionCC            Highs_metis_MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple         Highs_metis_MlevelNodeBisectionMultiple
#define MlevelNodeBisectionL2               Highs_metis_MlevelNodeBisectionL2
#define MlevelNodeBisectionL1               Highs_metis_MlevelNodeBisectionL1
#define SplitGraphOrder                     Highs_metis_SplitGraphOrder
#define SplitGraphOrderCC                   Highs_metis_SplitGraphOrderCC
#define MMDOrder                            Highs_metis_MMDOrder

/* options.c */
#define SetupCtrl                           Highs_metis_SetupCtrl
#define Setup2WayBalMultipliers             Highs_metis_Setup2WayBalMultipliers
#define PrintCtrl                           Highs_metis_PrintCtrl
#define CheckParams                         Highs_metis_CheckParams
#define FreeCtrl                            Highs_metis_FreeCtrl

/* refine.c */
#define Allocate2WayPartitionMemory         Highs_metis_Allocate2WayPartitionMemory
#define Compute2WayPartitionParams          Highs_metis_Compute2WayPartitionParams

/* separator.c */
#define ConstructSeparator                  Highs_metis_ConstructSeparator

/* sfm.c */
#define FM_2WayNodeRefine2Sided             Highs_metis_FM_2WayNodeRefine2Sided
#define FM_2WayNodeRefine1Sided             Highs_metis_FM_2WayNodeRefine1Sided
#define FM_2WayNodeBalance                  Highs_metis_FM_2WayNodeBalance

/* srefine.c */
#define Refine2WayNode                      Highs_metis_Refine2WayNode
#define Allocate2WayNodePartitionMemory     Highs_metis_Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams      Highs_metis_Compute2WayNodePartitionParams
#define Project2WayNodePartition            Highs_metis_Project2WayNodePartition

/* util.c */
#define iargmax_nrm                         Highs_metis_iargmax_nrm
#define iargmax2_nrm                        Highs_metis_iargmax2_nrm

/* wspace.c */
#define AllocateWorkSpace                   Highs_metis_AllocateWorkSpace
#define FreeWorkSpace                       Highs_metis_FreeWorkSpace
#define wspacemalloc                        Highs_metis_wspacemalloc
#define wspacepush                          Highs_metis_wspacepush
#define wspacepop                           Highs_metis_wspacepop
#define iwspacemalloc                       Highs_metis_iwspacemalloc
#define rwspacemalloc                       Highs_metis_rwspacemalloc
#define ikvwspacemalloc                     Highs_metis_ikvwspacemalloc


/* GKlib macros */
#define iaxpy                               Highs_metis_iaxpy
#define isum                                Highs_metis_isum
#define imalloc                             Highs_metis_imalloc
#define irealloc                            Highs_metis_irealloc
#define ismalloc                            Highs_metis_ismalloc
#define iset                                Highs_metis_iset
#define icopy                               Highs_metis_icopy
#define rmalloc                             Highs_metis_rmalloc
#define rrealloc                            Highs_metis_rrealloc
#define rsmalloc                            Highs_metis_rsmalloc
#define rset                                Highs_metis_rset
#define rcopy                               Highs_metis_rcopy
#define ikvmalloc                           Highs_metis_ikvmalloc
#define ikvrealloc                          Highs_metis_ikvrealloc
#define ikvsmalloc                          Highs_metis_ikvsmalloc
#define ikvset                              Highs_metis_ikvset
#define ikvcopy                             Highs_metis_ikvcopy
#define rkvmalloc                           Highs_metis_rkvmalloc
#define rkvrealloc                          Highs_metis_rkvrealloc
#define rkvsmalloc                          Highs_metis_rkvsmalloc
#define rkvset                              Highs_metis_rkvset
#define rkvcopy                             Highs_metis_rkvcopy
#define rpqCreate                           Highs_metis_rpqCreate
#define rpqInit                             Highs_metis_rpqInit
#define rpqReset                            Highs_metis_rpqReset
#define rpqFree                             Highs_metis_rpqFree
#define rpqDestroy                          Highs_metis_rpqDestroy
#define rpqLength                           Highs_metis_rpqLength
#define rpqInsert                           Highs_metis_rpqInsert
#define rpqDelete                           Highs_metis_rpqDelete
#define rpqUpdate                           Highs_metis_rpqUpdate
#define rpqGetTop                           Highs_metis_rpqGetTop
#define rpqSeeTopVal                        Highs_metis_rpqSeeTopVal
#define rpqSeeTopKey                        Highs_metis_rpqSeeTopKey
#define irand                               Highs_metis_irand
#define irandInRange                        Highs_metis_irandInRange
#define irandArrayPermute                   Highs_metis_irandArrayPermute
#define isortd                              Highs_metis_isortd
#define ikvsorti                            Highs_metis_ikvsorti

#endif