#ifndef HTIMERPRE_H_
#define HTIMERPRE_H_

#include <sys/time.h>
#include <cstdlib>
#include <cstdio>
#include <string>
using namespace std;

enum HTickItemPre {
   	// Presolve items
	HTICK_PRE_EMPTY_ROW,
	HTICK_PRE_FIXED,
	HTICK_PRE_SING_ROW,
	HTICK_PRE_DOUBLETON_EQUATION,
	HTICK_PRE_FORCING_ROW,
	HTICK_PRE_REDUNDANT_ROW,
	HTICK_PRE_DOMINATED_ROW_BOUNDS,
	HTICK_PRE_FREE_SING_COL,
	HTICK_PRE_SING_COL_DOUBLETON_INEQ,
	HTICK_PRE_IMPLIED_FREE_SING_COL,
	HTICK_PRE_DOMINATED_COLS,
	HTICK_PRE_WEAKLY_DOMINATED_COLS,
	HTICK_PRE_DOMINATED_COL_BOUNDS,
	HTICK_PRE_EMPTY_COL,
	//HTICK_PRE_DUPLICATE_ROWS,
	//HTICK_PRE_DUPLICATE_COLUMNS,
    // The total count
    HTICK_ITEMS_COUNT_PRE
};

class HTimerPre {
public:
    HTimerPre() {
    	itemNames.resize(HTICK_ITEMS_COUNT_PRE);

    	itemNames[HTICK_PRE_EMPTY_ROW] = "emptR";
    	itemNames[HTICK_PRE_FIXED] = "fixCo";
    	itemNames[HTICK_PRE_DOUBLETON_EQUATION] = "douEq";
        itemNames[HTICK_PRE_SING_ROW] = "singR";
        itemNames[HTICK_PRE_FORCING_ROW] = "forcR";
        itemNames[HTICK_PRE_REDUNDANT_ROW] = "reduR";
        itemNames[HTICK_PRE_FREE_SING_COL] = "colSi";
        itemNames[HTICK_PRE_SING_COL_DOUBLETON_INEQ] = "colSD";
        itemNames[HTICK_PRE_IMPLIED_FREE_SING_COL] = "colSI";
        itemNames[HTICK_PRE_DOMINATED_COLS] = "domCo";
        itemNames[HTICK_PRE_WEAKLY_DOMINATED_COLS] = "wdomC";
        itemNames[HTICK_PRE_EMPTY_COL] = "emptC";

        //itemNames[HTICK_PRE_DUPLICATE_ROWS] = "duplR";
        //itemNames[HTICK_PRE_DUPLICATE_COLUMNS] = "duplC";

        reset();
    }

    // Reset the timer
    void reset() {
        for (int i = 0; i < HTICK_ITEMS_COUNT_PRE; i++)
            itemTicks[i] = 0;
        startTime = getWallTime();
        startTick = getWallTick();
    }

    // Record tick of one item
    void recordStart(int itemKey) {
        itemStart[itemKey] = getWallTime();
    }
    void recordFinish(int itemKey) {
        itemTicks[itemKey] += (getWallTime() - itemStart[itemKey]);
    }

    // Report specific items in the list
    void report(int itemCount, int *itemList) {
        double totalTick = getTick();
        printf("Presolve times ");
        for (int i = 0; i < itemCount; i++)
            printf(" %s", itemNames[itemList[i]].c_str());
        printf("\n");
        printf("profile-item ");
        for (int i = 0; i < itemCount; i++) {
            int percent = 1000.0 * itemTicks[itemList[i]] / totalTick;
            printf(" %5d", percent);
        }
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

    double itemTicks[HTICK_ITEMS_COUNT_PRE];
    vector<string> itemNames;

private:
    // The start tick and time
    double startTime;
    double startTick;

    // The items tick and name
    double itemStart[HTICK_ITEMS_COUNT_PRE];


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

#endif /* HTIMERPRE_H_ */
