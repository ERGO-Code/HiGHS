// forward declare a struct called MyThing, but make sure the contents of it aren't available so that we can do anything we want to it in the implementation.
struct MyThing;

// define a class that contains function pointers.  We have a requirement that these function pointers must be valid for the lifetime of them being registered with MyThing.
struct MyThingCallback {
    void (*eventA)(MyThing* sender, const char* someData);
    void (*eventB)(MyThing* sender, int someOtherData);
};

// Some helper functions to do things with MyThing
MyThing* createMyThing();
void deleteMyThing(MyThing*);
void registerCallback(MyThing*, MyThingCallback);

int getSomeProperty(const MyThing*);
void setSomeProperty(MyThing*, int);

class Highs;
Highs* createHighs();
void deleteHighs(Highs*);