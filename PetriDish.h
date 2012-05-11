//
//  PetriDish.h
//  PetriDish v0.9
//
//  Replaced by Glenn Sugden on 2012.05.10.
//  Created by Glenn Sugden on 2012.05.05.
//  This source is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//  To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
//

#ifndef GALib_PetriDish_h
#define GALib_PetriDish_h

//================================================================================
#pragma mark Definitions
//================================================================================

#define GA_CLASS            GASimpleGA     // GASimpleGA     // GADemeGA  // GASteadyStateGA // GAIncrementalGA
#define USING_GASIMPLEGA    1              // Set this to 0 if using one of the other GA's above

//================================================================================
#pragma mark Headers
//================================================================================

#include <ga/ga.h>

//================================================================================
#pragma mark Classes
//================================================================================

class PetriDish : public GA_CLASS
{
    
public:
    
    // CONSTRUCTOR
    PetriDish(const GAGenome& genomeToClone);
    
    // EVALUATOR
    static void GAPDEvaluator( GAPopulation& pop );
    
    // TERMINATOR (interrupt)
    static GABoolean InterruptTerminator(GAGeneticAlgorithm & ga);

    // DESTRUCTOR
    ~PetriDish();
    
    // UTILITY
    void interrupt()     { _interrupt = true; }
    
private:
    
    bool _interrupt;            // A flag indicating we need to bail on the Evaluation
    
#if USING_GASIMPLEGA
    bool _oldPopInitialized;    // A flag for a GASimpleGA workaround (see GAPDEvaluator)
#endif//USING_GASIMPLEGA
    
};

#endif
