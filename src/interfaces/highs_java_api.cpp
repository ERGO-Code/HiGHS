//
// Created by poteh on 22.04.2020.
//
#include "Highs.h"
#include "highs_java_api.h"
#include "highs_c_api.h"
#include "stdio.h"

#include <chrono>
#include <thread>


// Пример экспортированной функции.
JNIEXPORT jint JNICALL Java_org_best_team_Highs_multiply
        (JNIEnv *env, jobject, jint a, jint b) {
return a*b;
}

void printIntArray(const string& varName, int* array, int size);

void printDoubleArray(const string& varName, double* array, int size);

// Вызов процедуры линейной оптимизации.
JNIEXPORT jobject JNICALL Java_org_best_team_Highs_invokeLpOptimization
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

    // Read java reflection metadata
    jclass thisObjClass = env->GetObjectClass(thisObj);
    jmethodID setStatusMethod = env->GetMethodID(thisObjClass,"setStatus", "(I)V");
    jmethodID onCompleteMethod = env->GetMethodID(thisObjClass,"onComplete", "(I)V");

    jobject solutionObject = nullptr;
    if(env->ExceptionCheck()){
        std::cout << "Exception occurred" << std::endl;
        return solutionObject;
    }


    jclass solutionClass = env->FindClass("org/best/team/HighsSolution");
    jmethodID constructorId = env->GetMethodID(solutionClass, "<init>", "(II)V");
    solutionObject = env->NewObject(solutionClass, constructorId, 10, 10);


    //   jfieldID solutionFieldId = env->GetFieldID(modelObjClass,
    //           "solution", "Lorg/best/team/HighsSolution;");

    //   if(((long)solutionFieldId) == 0){
    //      std::cout << "solution id pointer is not set" << std::endl;
    //   } else {
//        std::cout << "solution field ID was found" << std::endl;
//    }

 //   env -> SetObjectField(modelObj, solutionFieldId, solutionObject);



    jdouble *colcost = env->GetDoubleArrayElements(jColcost, JNI_FALSE);
    jdouble *collower = env->GetDoubleArrayElements(jCollower, JNI_FALSE);
    jdouble *colupper = env->GetDoubleArrayElements(jColupper, JNI_FALSE);
    jdouble *rowlower = env->GetDoubleArrayElements(jRowlower, JNI_FALSE);
    jdouble *rowupper = env->GetDoubleArrayElements(jRowupper, JNI_FALSE);
    jint *astart = env->GetIntArrayElements(jAstart, JNI_FALSE);
    jint *aindex = env->GetIntArrayElements(jAindex, JNI_FALSE);
    jdouble *avalue = env->GetDoubleArrayElements(jAvalue, JNI_FALSE);


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


    if(env->ExceptionCheck()){
        std::cout << "Exception occurred" << std::endl;
        return solutionObject;
    }

    Highs highs;

    int status =
            Highs_passLp(&highs, numcol, numrow, numnz, colcost, collower, colupper,
                        rowlower, rowupper, astart, aindex, avalue);

    std::cout << "Called Highs_passLp with status " << status << std::endl;

    env -> CallVoidMethod(thisObj, setStatusMethod, status);
  //  std::this_thread::sleep_for(std::chrono::nanoseconds(1000000));


    std::cout << "Trying to call highs.run" << std::endl;

    status = (int)highs.run();


    std::cout << "highs.run with status" << status << std::endl;

    if (status == 0) {
     //   HighsSolution solution;
    //    HighsBasis basis;
    //    solution = highs.getSolution();
     //   basis = highs.getBasis();
        int modelStatus = (int)highs.getModelStatus();
         std::cout << "model status is "<< modelStatus << std::endl;

        for (int i = 0; i < numcol; i++) {

           // std::cout << "column "<< i << "=" << solution.col_value[i] << std::endl;

          //  colvalue[i] = solution.col_value[i];
         //   coldual[i] = solution.col_dual[i];

         //   colbasisstatus[i] = (int)basis.col_status[i];
        }

        for (int i = 0; i < numrow; i++) {
         //   rowvalue[i] = solution.row_value[i];
        //    rowdual[i] = solution.row_dual[i];

          //  rowbasisstatus[i] = (int)basis.row_status[i];
        }
    } else {
        std::cout << "error" << std::endl;
    }


  //  env -> CallVoidMethod(thisObj, setStatusMethod, status);


   // std::this_thread::sleep_for(std::chrono::nanoseconds(1000000));


   // env -> CallVoidMethod(thisObj, setStatusMethod, status);


    env -> CallVoidMethod(thisObj, onCompleteMethod, status);

    env->ReleaseDoubleArrayElements(jColcost, colcost,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jCollower, collower,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jCollower, colupper,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jRowlower, rowlower,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jRowupper, rowupper,  JNI_ABORT);
    env->ReleaseIntArrayElements(jAstart, astart,  JNI_ABORT);
    env->ReleaseIntArrayElements(jAindex, aindex,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(jAvalue, avalue,  JNI_ABORT);

    return solutionObject;

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