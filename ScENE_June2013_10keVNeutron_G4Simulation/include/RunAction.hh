#ifndef __RUNACTION_H__
#define __RUNACTION_H__

#include <G4UserRunAction.hh>
#include "AnalysisManager.hh"


class G4Run;
class AnalysisManager;


class RunAction: public G4UserRunAction
{
    public:
    RunAction(AnalysisManager *AnalysisManagerPointer){
        pAnalysisManager = AnalysisManagerPointer;
    }
        ~RunAction();

    public:
    void BeginOfRunAction(const G4Run *pRun);
    void   EndOfRunAction(const G4Run *pRun);
    
    private:
    AnalysisManager *pAnalysisManager;
};

#endif __RUNACTION_H__