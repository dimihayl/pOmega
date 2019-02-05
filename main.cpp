
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "pOmega_1.h"

#include "DLM_CppTools.h"

using namespace std;


int main(int argc, char *argv[])
{

	printf("Executing the pOmega main...\n");
    DLM_Timer TIMER;

    SimpleFitter1();

    char** ARGV = NULL;
    if(argc) ARGV=new char* [argc];
    for(int iARG=1; iARG<argc; iARG++){
        ARGV[iARG] = new char [128];
        strcpy(ARGV[iARG],argv[iARG]);
    }

    for(int iARG=1; iARG<argc; iARG++){
        delete [] ARGV[iARG];
    }
    if(ARGV) delete [] ARGV;

    long long ExeTime = TIMER.Stop()/1000.;
    char* strtime = new char [128];
    ShowTime(ExeTime,strtime,0,true,6);
    printf("The script terminated after: %s\n",strtime);

    delete [] strtime;


    return 0;
}
