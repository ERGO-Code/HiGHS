#ifndef HTIMER_H_
#define HTIMER_H_

#include <sys/time.h>
#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <string>
using namespace std;

enum HTickItem {
    // Predefined groups
    HTICK_GROUP1 = 0,
    HTICK_GROUP2,
    HTICK_GROUP3,
    HTICK_GROUP4,
    HTICK_GROUP5,
    HTICK_GROUP6,
    HTICK_GROUP7,
    HTICK_GROUP8,
    HTICK_GROUP9,

    // Simplex method related items
    HTICK_INVERT,
    HTICK_CHUZR1,
    HTICK_CHUZC0,
    HTICK_CHUZC1,
    HTICK_CHUZC2,
    HTICK_CHUZC3,
    HTICK_CHUZC4,
    HTICK_FTRAN,
    HTICK_BTRAN,
    HTICK_PRICE,
    HTICK_FTRAN_DSE,
    HTICK_FTRAN_MIX,
    HTICK_UPDATE_DUAL,
    HTICK_UPDATE_PRIMAL,
    HTICK_UPDATE_WEIGHT,
    HTICK_UPDATE_FACTOR,
    HTICK_UPDATE_ROW_EP,

    // The total count
    HTICK_ITEMS_COUNT
};

class HTimer {
public:
  HTimer() {
    itemNames[HTICK_GROUP1] = "GROUP1"; itemCh3Names[HTICK_GROUP1] = "GP1";
    itemNames[HTICK_GROUP2] = "GROUP2"; itemCh3Names[HTICK_GROUP2] = "GP2";
    itemNames[HTICK_GROUP3] = "GROUP3"; itemCh3Names[HTICK_GROUP3] = "GP3";
    itemNames[HTICK_GROUP4] = "GROUP4"; itemCh3Names[HTICK_GROUP4] = "GP4";
    itemNames[HTICK_GROUP5] = "GROUP5"; itemCh3Names[HTICK_GROUP5] = "GP5";
    itemNames[HTICK_GROUP6] = "GROUP6"; itemCh3Names[HTICK_GROUP6] = "GP6";
    itemNames[HTICK_GROUP7] = "GROUP7"; itemCh3Names[HTICK_GROUP7] = "GP7";
    itemNames[HTICK_GROUP8] = "GROUP8"; itemCh3Names[HTICK_GROUP8] = "GP8";
    itemNames[HTICK_GROUP9] = "GROUP9"; itemCh3Names[HTICK_GROUP9] = "GP9";

    itemNames[HTICK_INVERT] = "INVERT"; itemCh3Names[HTICK_INVERT] = "INV";
    itemNames[HTICK_CHUZR1] = "CHUZR1"; itemCh3Names[HTICK_CHUZR1] = "CR1";
    itemNames[HTICK_CHUZC0] = "CHUZC0"; itemCh3Names[HTICK_CHUZC0] = "CC0";
    itemNames[HTICK_CHUZC1] = "CHUZC1"; itemCh3Names[HTICK_CHUZC1] = "CC1";
    itemNames[HTICK_CHUZC2] = "CHUZC2"; itemCh3Names[HTICK_CHUZC2] = "CC2";
    itemNames[HTICK_CHUZC3] = "CHUZC3"; itemCh3Names[HTICK_CHUZC3] = "CC3";
    itemNames[HTICK_CHUZC4] = "CHUZC4"; itemCh3Names[HTICK_CHUZC4] = "CC4";
    itemNames[HTICK_FTRAN] = "FTRAN"; itemCh3Names[HTICK_FTRAN] = "COL";
    itemNames[HTICK_BTRAN] = "BTRAN"; itemCh3Names[HTICK_BTRAN] = "REP";
    itemNames[HTICK_PRICE] = "PRICE"; itemCh3Names[HTICK_PRICE] = "RAP";
    itemNames[HTICK_FTRAN_DSE] = "FTRAN_DSE"; itemCh3Names[HTICK_FTRAN_DSE] = "DSE";
    itemNames[HTICK_FTRAN_MIX] = "FTRAN_MIX"; itemCh3Names[HTICK_FTRAN_MIX] = "MIX";
    itemNames[HTICK_UPDATE_DUAL] = "UPDATE_DUAL"; itemCh3Names[HTICK_UPDATE_DUAL] = "UPD";
    itemNames[HTICK_UPDATE_PRIMAL] = "UPDATE_PRIMAL"; itemCh3Names[HTICK_UPDATE_PRIMAL] = "UPP";
    itemNames[HTICK_UPDATE_WEIGHT] = "UPDATE_WEIGHT"; itemCh3Names[HTICK_UPDATE_WEIGHT] = "UPW";
    itemNames[HTICK_UPDATE_FACTOR] = "UPDATE_FACTOR"; itemCh3Names[HTICK_UPDATE_FACTOR] = "UPF";
    itemNames[HTICK_UPDATE_ROW_EP] = "UPDATE_ROW_EP"; itemCh3Names[HTICK_UPDATE_ROW_EP] = "UPR";
    
    reset();
  }

  // Reset the timer
  void reset() {
    for (int i = 0; i < HTICK_ITEMS_COUNT; i++) {
      itemNumCall[i] = 0;
      itemStart[i] = initialItemStart;
      itemTicks[i] = 0;
    }
    startTime = getWallTime();
    startTick = getWallTick();
  }

  // Record tick of one item
  void recordStart(int itemKey) {
    // Check that the clock's been stopped: should be set to
    // getWallTick() >= 0 [or initialised to initialItemStart > 0]
    //    if (itemKey == HTICK_CHUZR1 || itemKey == HTICK_INVERT)
    if (itemStart[itemKey]<0)
      printf("recordStart [%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
	     itemKey, itemStart[itemKey], itemTicks[itemKey], itemNumCall[itemKey]);
    assert(itemStart[itemKey]<0);
    // Set the start to be the negation of the WallTick to check that
    // the clock's been started in recordFinish
    itemStart[itemKey] = -getWallTick();
  }
  void recordFinish(int itemKey) {
    // Check that the clock's been started: should be set to
    // -getWallTick() <= 0
    //    if (itemKey == HTICK_CHUZR1 || itemKey == HTICK_INVERT)
    if (itemStart[itemKey]>0)
      printf("recordFinish[%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
	     itemKey, itemStart[itemKey], itemTicks[itemKey], itemNumCall[itemKey]);
    assert(itemStart[itemKey]>0);
    double wallTick = getWallTick();
    itemTicks[itemKey] += (wallTick + itemStart[itemKey]);
    itemNumCall[itemKey]++;
    // Set the start to be the WallTick to check that the clock's been
    // stopped in recordStart
    itemStart[itemKey] = wallTick;
  }

  // Report specific items in the list
  void report(int itemCount, int *itemList) {
    double totalTick = getTick();
    printf("txt-profile-name ");
    for (int i = 0; i < itemCount; i++)
      printf(" %-3s", itemCh3Names[itemList[i]].c_str());
    printf("\n");
    printf("txt-profile-item ");
    int suPerMille=0;
    for (int i = 0; i < itemCount; i++) {
      int perMille = 1000.0 * itemTicks[itemList[i]] / totalTick;
      printf(" %3d", perMille);
      suPerMille += perMille;
    }
    printf(" per mille: Sum = %d", suPerMille);
    printf("\n");
    printf("txt-profile-time ID: Operation       :     Time           :    Calls  Time/Call\n");
    double tick2sec = 3.6e-10;
    double suTi=0;
    for (int i = 0; i < itemCount; i++) {
      int item = itemList[i];
      double ti = tick2sec*itemTicks[item];
      double perCent = 100.0 * itemTicks[item] / totalTick;
      double tiPerCall = 0;
      if (itemNumCall[item]>0) tiPerCall = ti/itemNumCall[item];
      printf("txt-profile-time %2d: %-16s: %11.4e (%4.1f%%): %7d %11.4e\n",
	     item, itemNames[item].c_str(), ti, perCent, itemNumCall[item], tiPerCall);
      suTi += ti;
    }
    double perCent = 100.0 * suTi / totalTick;
    printf("txt-profile-time   : SUM             : %11.4e (%4.1f%%)\n", suTi, perCent);
    printf("txt-profile-time   : TOTAL           : %11.4e\n", tick2sec*totalTick);
    //Report for Excel
    printf("grep_excel-profile-name");
    for (int i = 0; i < itemCount; i++) printf(",%s", itemNames[itemList[i]].c_str());
    printf(",SumTime,");
    printf(",TotalTime\n");
    printf("grep_excel-profile-time");
    for (int i = 0; i < itemCount; i++) printf(",%e", tick2sec*itemTicks[itemList[i]]);
    printf(",%e,", suTi);
    printf(",%e\n", tick2sec*totalTick);
    printf("grep_excel-profile-calls");
    for (int i = 0; i < itemCount; i++) printf(",%d", itemNumCall[itemList[i]]);
    printf("\n");
  }
  
  // Current eclipsed time
  double getTime() {
    return getWallTime() - startTime;
  }
  
  // Current eclipsed CPU tick
  double getTick() {
    return getWallTick() - startTick;
  }
  
  // private:
  // The start tick and time
  double startTime;
  double startTick;

  // The items tick and name
  const double initialItemStart = 1.0;
  int itemNumCall[HTICK_ITEMS_COUNT];
  double itemStart[HTICK_ITEMS_COUNT];
  double itemTicks[HTICK_ITEMS_COUNT];
  string itemNames[HTICK_ITEMS_COUNT];
  string itemCh3Names[HTICK_ITEMS_COUNT];
  
  // Current wall time
  double getWallTime() {
    double walltime;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    walltime = tv.tv_sec;
    walltime += (double) tv.tv_usec / 1000000.0;
    return walltime;
  }
  
  // Current wall tick
  double getWallTick() {
    unsigned a, d;
    asm volatile("rdtsc" : "=a" (a), "=d" (d));
    return ((unsigned long long) a) | (((unsigned long long) d) << 32);
  }
};

#endif /* HTIMER_H_ */
