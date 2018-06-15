#ifndef HTIMER_H_
#define HTIMER_H_

#include <sys/time.h>
#include <cstdlib>
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
        itemNames[HTICK_GROUP1] = "GP1";
        itemNames[HTICK_GROUP2] = "GP2";
        itemNames[HTICK_GROUP3] = "GP3";
        itemNames[HTICK_GROUP4] = "GP4";
        itemNames[HTICK_GROUP5] = "GP5";
        itemNames[HTICK_GROUP6] = "GP6";
        itemNames[HTICK_GROUP7] = "GP7";
        itemNames[HTICK_GROUP8] = "GP8";
        itemNames[HTICK_GROUP9] = "GP9";

        itemNames[HTICK_INVERT] = "INV";
        itemNames[HTICK_CHUZR1] = "CR1";
        itemNames[HTICK_CHUZC0] = "CC0";
        itemNames[HTICK_CHUZC1] = "CC1";
        itemNames[HTICK_CHUZC2] = "CC2";
        itemNames[HTICK_CHUZC3] = "CC3";
        itemNames[HTICK_CHUZC4] = "CC4";
        itemNames[HTICK_FTRAN] = "COL";
        itemNames[HTICK_BTRAN] = "REP";
        itemNames[HTICK_PRICE] = "RAP";
        itemNames[HTICK_FTRAN_DSE] = "DSE";
        itemNames[HTICK_FTRAN_MIX] = "MIX";
        itemNames[HTICK_UPDATE_DUAL] = "UPD";
        itemNames[HTICK_UPDATE_PRIMAL] = "UPP";
        itemNames[HTICK_UPDATE_WEIGHT] = "UPW";
        itemNames[HTICK_UPDATE_FACTOR] = "UPF";
        itemNames[HTICK_UPDATE_ROW_EP] = "UPR";

        reset();
    }

    // Reset the timer
    void reset() {
        for (int i = 0; i < HTICK_ITEMS_COUNT; i++)
            itemTicks[i] = 0;
        startTime = getWallTime();
        startTick = getWallTick();
    }

    // Record tick of one item
    void recordStart(int itemKey) {
        itemStart[itemKey] = getWallTick();
    }
    void recordFinish(int itemKey) {
        itemTicks[itemKey] += (getWallTick() - itemStart[itemKey]);
    }

    // Report specific items in the list
    void report(int itemCount, int *itemList) {
        double totalTick = getTick();
        printf("profile-name ");
        for (int i = 0; i < itemCount; i++)
            printf(" %s", itemNames[itemList[i]].c_str());
        printf("\n");
        printf("profile-item ");
	int suPerMille=0;
        for (int i = 0; i < itemCount; i++) {
            int perMille = 1000.0 * itemTicks[itemList[i]] / totalTick;
            printf(" %3d", perMille);
	    suPerMille += perMille;
        }
	printf(" per mille: Sum = %d", suPerMille);
        printf("\n");
	//	        printf("profile-time\n");
	//		double tick2sec = 3.6e-10;
	//		double suTi=0;
	//	        for (int i = 0; i < itemCount; i++) {
	//		  double ti = tick2sec*itemTicks[itemList[i]];
	//		  printf(" %s: %6.4g\n", itemNames[itemList[i]].c_str(), ti);
	//		  suTi += ti;
	//	        }
	//		printf(" suTi = %g\n", suTi);
	//		printf(" totalTick = %g\n", tick2sec*totalTick);
    }

    void reportLocal(int itemCount, int *itemList) {
        double totalTick = 0;
        for (int i = 0; i < itemCount; i++) {
            totalTick += itemTicks[itemList[i]];
        }
        printf("profile-name ");
        for (int i = 0; i < itemCount; i++)
            printf(" %s", itemNames[itemList[i]].c_str());
        printf("\n");
        printf("profile-item ");
	int suPerMille=0;
        for (int i = 0; i < itemCount; i++) {
            int perMille = 1000.0 * itemTicks[itemList[i]] / totalTick;
            printf(" %3d", perMille);
	    suPerMille += perMille;
        }
	printf(" per mille: Sum = %d", suPerMille);
        printf("\n");

	//        printf("profile-time\n");
	//        for (int i = 0; i < itemCount; i++) {
	//	  printf(" %s: %6.4g\n", itemNames[itemList[i]].c_str(), itemTicks[itemList[i]]);
	//        }
	//	printf(" totalTick = %g\n", totalTick);

    }

    // Current eclipsed time
    double getTime() {
        return getWallTime() - startTime;
    }

    // Current eclipsed CPU tick
    double getTick() {
        return getWallTick() - startTick;
    }

private:
    // The start tick and time
    double startTime;
    double startTick;

    // The items tick and name
    double itemStart[HTICK_ITEMS_COUNT];
    double itemTicks[HTICK_ITEMS_COUNT];
    string itemNames[HTICK_ITEMS_COUNT];

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
