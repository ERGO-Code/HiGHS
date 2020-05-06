//
// Created by potehin on 22.04.2020.
//
#include "Highs.h"
#include "highs_java_api.h"
#include "highs_c_api.h"

#include <thread>

void printIntArray(const string& varName, int* array, int size);

void printDoubleVector(const string& varName,  vector<double> vector);

void printStatusVector(const string& varName, vector<HighsBasisStatus> vector);

void printDoubleArray(const string& varName, double* array, int size);

// Вызов процедуры линейной оптимизации.
JNIEXPORT void JNICALL Java_org_best_team_Highs_invokeLpOptimization
        (JNIEnv *env, jobject thisObj,
         jint numcol,        //!< number of columns
         jint numrow,        //!< number of rows
         jint numnz,       //!< number of entries in the constraint matrix
         jdoubleArray jColcost,   //!< array of length [numcol] with column costs
         jdoubleArray jCollower,  //!< array of length [numcol] with lower column bounds
         jdoubleArray jColupper,  //!< array of length [numcol] with upper column bounds
         jdoubleArray jRowlower,  //!< array of length [numrow] with lower row bounds
         jdoubleArray jRowupper,  //!< array of length [numrow] with upper row bounds
         jintArray jAstart,       //!< array of length [numcol+1] with column start indices
         jintArray jAindex,  //!< array of length [numnz] with row indices of matrix entries
         jdoubleArray jAvalue  //!< array of length [numnz] with value of matrix entries
         ) {

    std::cout << "Call LpOptimization from cpp" << std::endl;

    bool toPrintModel = numcol <10 && numrow <10; // print only small models for debug

    // Read java reflection metadata
    jclass thisObjClass = env->GetObjectClass(thisObj);
    jmethodID setStatusMethod = env->GetMethodID(thisObjClass,"setStatus", "(I)V");
    jmethodID onCompleteMethod = env->GetMethodID(thisObjClass,
            "onComplete", "(IILorg/best/team/HighsSolution;Lorg/best/team/HighsBasis;)V");


    if(env->ExceptionCheck()){  //check if deferred  exception occurred
        std::cout << "Exception occurred" << std::endl;
        return;
    }


    jclass solutionClass = env->FindClass("org/best/team/HighsSolution");
    jmethodID solutionConstructorId = env->GetMethodID(solutionClass, "<init>", "(II)V");
    jobject solutionObject = env->NewObject(solutionClass, solutionConstructorId, numcol, numrow);

    auto colValueArray = (jdoubleArray)env->GetObjectField(solutionObject,
                                                      env->GetFieldID(solutionClass,"colvalue", "[D")
                                                      );
    auto colDualArray = (jdoubleArray)env->GetObjectField(solutionObject,
                                                      env->GetFieldID(solutionClass,"coldual", "[D")
    );
    auto rowValueArray = (jdoubleArray)env->GetObjectField(solutionObject,
                                                      env->GetFieldID(solutionClass,"rowvalue", "[D")
    );
    auto rowDualArray = (jdoubleArray)env->GetObjectField(solutionObject,
                                                      env->GetFieldID(solutionClass,"rowdual", "[D")
    );


    jclass basisClass = env->FindClass("org/best/team/HighsBasis");
    jmethodID basisConstructorId = env->GetMethodID(basisClass, "<init>", "(II)V");
    jobject basisObject = env->NewObject(basisClass, basisConstructorId, numcol, numrow);

    auto colBasisStatusArray = (jintArray)env->GetObjectField(basisObject,
                                                      env->GetFieldID(basisClass,"colbasisstatus", "[I")
    );
    auto rowBasisStatusArray = (jintArray)env->GetObjectField(basisObject,
                                                     env->GetFieldID(basisClass,"rowbasisstatus", "[I")
    );

    jdouble *colcost = env->GetDoubleArrayElements(jColcost, JNI_FALSE);
    jdouble *collower = env->GetDoubleArrayElements(jCollower, JNI_FALSE);
    jdouble *colupper = env->GetDoubleArrayElements(jColupper, JNI_FALSE);
    jdouble *rowlower = env->GetDoubleArrayElements(jRowlower, JNI_FALSE);
    jdouble *rowupper = env->GetDoubleArrayElements(jRowupper, JNI_FALSE);
    jint *astart = env->GetIntArrayElements(jAstart, JNI_FALSE);
    jint *aindex = env->GetIntArrayElements(jAindex, JNI_FALSE);
    jdouble *avalue = env->GetDoubleArrayElements(jAvalue, JNI_FALSE);

    if(env->ExceptionCheck()){ //check if deferred  exception occurred
        std::cout << "Exception occurred" << std::endl;
        return;
    }



    if(toPrintModel){  // print only small models
        std::cout << "numcol=" << numcol <<  std::endl;
        std::cout << "numrow=" << numrow <<  std::endl;
        std::cout << "numnz=" << numnz <<  std::endl;
        printDoubleArray("colcost", colcost, env->GetArrayLength(jColcost));
        printDoubleArray("collower", collower, env->GetArrayLength(jCollower));
        printDoubleArray("colupper", colupper, env->GetArrayLength(jColupper));
        printDoubleArray("rowlower", rowlower, env->GetArrayLength(jRowlower));
        printDoubleArray("rowupper", rowupper, env->GetArrayLength(jRowupper));
        printIntArray("astart", astart, env->GetArrayLength(jAstart));
        printIntArray("aindex", aindex, env->GetArrayLength(jAindex));
        printDoubleArray("avalue", avalue, env->GetArrayLength(jAvalue));
    }





    Highs highs;

    int status =
            Highs_passLp(&highs, numcol, numrow, numnz, colcost, collower, colupper,
                        rowlower, rowupper, astart, aindex, avalue);
    env -> CallVoidMethod(thisObj, setStatusMethod, status);

    status = (int)highs.run();

    if (status == 0) {
        HighsSolution solution;
        HighsBasis basis;
        solution = highs.getSolution();
        basis = highs.getBasis();
        int modelStatus = (int)highs.getModelStatus();

        if(toPrintModel) {
            std::cout << "model status is " << modelStatus << std::endl;
            printDoubleVector("colvalue", solution.col_value);
            printDoubleVector("coldual", solution.col_dual);
            printDoubleVector("rowvalue", solution.row_value);
            printDoubleVector("rowdual", solution.row_dual);
            printStatusVector("colbasisstatus", basis.col_status);
            printStatusVector("rowbasisstatus", basis.row_status);
        }

        auto colvalue = (double*)malloc(sizeof(double) * numcol);
        auto coldual = (double*)malloc(sizeof(double) * numcol);
        auto colbasisstatus = (int*)malloc(sizeof(int) * numcol);

        for (int i = 0; i < numcol; i++) {
            colvalue[i] = solution.col_value[i];
            coldual[i] = solution.col_dual[i];
            colbasisstatus[i] = (int)basis.col_status[i];
        }
        env -> SetDoubleArrayRegion(colValueArray, 0, solution.col_value.size(),colvalue);
        env -> SetDoubleArrayRegion(colDualArray, 0, solution.col_dual.size(),coldual);
        env -> SetIntArrayRegion(colBasisStatusArray, 0, basis.col_status.size(),colbasisstatus);

        free(colvalue);
        free(coldual);
        free(colbasisstatus);

        auto rowvalue = (double*)malloc(sizeof(double) * numrow);
        auto rowdual = (double*)malloc(sizeof(double) * numrow);
        auto rowbasisstatus = (int*)malloc(sizeof(int) * numrow);

        for (int i = 0; i < numrow; i++) {
            rowvalue[i] = solution.row_value[i];
            rowdual[i] = solution.row_dual[i];
            rowbasisstatus[i] = (int)basis.row_status[i];
        }
        env -> SetDoubleArrayRegion(rowValueArray, 0, solution.row_value.size(),rowvalue);
        env -> SetDoubleArrayRegion(rowDualArray, 0, solution.row_dual.size(),rowdual);
        env -> SetIntArrayRegion(rowBasisStatusArray, 0, basis.row_status.size(),rowbasisstatus);

        free(rowvalue);
        free(rowdual);
        free(rowbasisstatus);

        env -> CallVoidMethod(thisObj, onCompleteMethod, status, modelStatus, solutionObject, basisObject);
    } else {
        std::cout << "error" << std::endl;
    }

    env->ReleaseDoubleArrayElements(jColcost, colcost,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jCollower, collower,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jCollower, colupper,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jRowlower, rowlower,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jRowupper, rowupper,  JNI_ABORT);
    env->ReleaseIntArrayElements(jAstart, astart,  JNI_ABORT);
    env->ReleaseIntArrayElements(jAindex, aindex,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jAvalue, avalue,  JNI_ABORT);

}



void printDoubleVector(const string& varName, vector<double> vector){
    std::cout << varName << "=[";
    for (unsigned i=0; i<vector.size(); i++){
        std::cout << vector[i] ;
        if(i < vector.size() - 1){
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}


void printStatusVector(const string& varName, vector<HighsBasisStatus> vector){
    std::cout << varName << "=[";
    for (unsigned i=0; i<vector.size(); i++){
        std::cout << (int)vector[i];
        if(i < vector.size() - 1){
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

void printDoubleArray(const string& varName,  double* array, int size){
    std::cout << varName << "=[";
    for (int i=0; i<size; i++){
        std::cout << array[i] ;
        if(i < size - 1){
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

void printIntArray(const string& varName,  int* array, int size){
    std::cout << varName << "=[";
    for (int i=0; i<size; i++){
        std::cout << array[i] ;
        if(i < size - 1){
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}