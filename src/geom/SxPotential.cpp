// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxlib.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#include <SxPotential.h>
#include <SxGrid.h>

SxPotential::SxPotential ()
   : dEnergyLow(false), noNetForce(true)
{
   // empty
}

SxPotential::~SxPotential ()
{
   // empty
}

SxAtomicStructure
SxPotential::getForces (const SxAtomicStructure &tau,
                        const SxArray<const SxSymbolTable *> &cmds)
{
   SxAtomicStructure f;
   for (int i=0; i < cmds.getSize(); i++)
         f = getForces (tau, cmds(i));

   return f;
}

SxAtomicStructure SxPotential::getSymForces (const SxAtomicStructure  &tau,
                                             const SxSymbolTable      *table)
{
   SX_CHECK (isSymmetrizedStructure (tau));

   // --- compute unsymmetrized forces
   SxAtomicStructure f = getForces (tau, table);

   // symmetrize forces according to existing symmetry elements
   if (forceSymmetrizer.getNSymmetries () == 0) {
      forceSymmetrizer.setup (tau);
   }

   return forceSymmetrizer | f;
}

SxAtomicStructure
SxPotential::getSymForces (const SxAtomicStructure &tau,
                           const SxArray<const SxSymbolTable *> &cmds)
{
   if (cmds.getSize () == 0) {
      SX_CHECK (isRegistered (NULL));
      return getSymForces (tau, NULL);
   }
   SxAtomicStructure f;
   for (int i=0; i < cmds.getSize(); i++) {
      f = getSymForces (tau, cmds(i));
   }
   return f;
}


SxArray<const SxSymbolTable *>
SxPotential::getMinimCmds (const SxSymbolTable *table) const
{
   SX_CHECK (table);
   SxList<const SxSymbolTable *> cmdList;
   try  {
      if (table->containsGroup ("bornOppenheimer"))
         table = table->getGroup("bornOppenheimer");

      SxSymbolTable *cmd   = NULL;
      for (cmd  = table->begin();
           cmd != NULL;
           cmd  = cmd->nextSibling())
      {
         if (isRegistered (cmd))  cmdList << cmd;
         else  {
            if (!table->useDeprecateFormat())  {
               sxprintf ("Unknown command.\n");
               cmd->print ();
               SX_EXIT;
            }
         }
      }

   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
   return SxArray<const SxSymbolTable *> (cmdList);
}


bool SxPotential::isSymmetrizedStructure (const SxAtomicStructure &tau) const
{
   SX_CHECK (tau.cell.symGroupPtr);
   const SxSymGroup &S = *tau.cell.symGroupPtr;

   int nSym = S.getSize();

   // --- verify existing symmetry elements
   if (tau.getNAtoms () < 100)  {
      for (int iSym=0; iSym < nSym; iSym++)
         if (tau != (S(iSym) ^ tau) )
            return false;
   } else {
      SxGrid grid(tau, 3);
      for (int iSym=0; iSym < nSym; iSym++)
         if (! tau.match (grid, S(iSym) ^ tau) )
            return false;
   }

   if (forceSymmetrizer.getNSymmetries () != 0)
      if (!forceSymmetrizer.checkStr (tau)) return false;
   return true;

}

