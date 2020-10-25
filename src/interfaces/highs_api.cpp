#include <interfaces/highs_api.h>
#include <src/Highs.h>

#include <vector>

struct MyThing {
private:
    int property;
    std::vector<MyThingCallback> callbacks;
    // functions are allowed here
};

Highs* createHighs() {
  Highs* highs = new Highs();
}

void deleteHighs(Highs*);