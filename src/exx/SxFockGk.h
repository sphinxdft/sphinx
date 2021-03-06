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

#ifndef _SX_FOCK_GK_H_
#define _SX_FOCK_GK_H_

#include <SxExx.h>
#include <SxGBasis.h>
#include <SxGkBasis.h>
#include <SxRBasis.h>
#include <SxDirac.h>
#include <SxTypes.h>
#include <SxCell.h>
#include <SxPW.h>
#include <SxFermi.h>
#include <SxArray.h>

/** \brief Definition of the Fock operator used in EXX calculations (old,
           non-symmetrised version, see SxSymFockGk for the symmetrised
           version)

    \b SxFockGkk = S/PHI/nX Fock Operator in reciprocal space

    Computes the exchange (Fock) operator \f$ V_{\mathbf k}^{\rm x} \f$
    in reciprocal space:
    \f[
      V_{\rm x}({\mathbf k}, {\mathbf G}, {\mathbf G'})
        = -\frac{4\pi {\rm e}^2}{\Omega_{\rm cell}}
          \sum_{v, {\mathbf q}, {\mathbf G}_1}
          \frac{c_{v\mathbf q}({\mathbf G}+{\mathbf G}_1)
                c_{v\mathbf q}^{\ast}({\mathbf G'}+{\mathbf G}_1)}
               {|{\mathbf q} - {\mathbf k} + {\mathbf G}_1|^2}
    \f]

    \f$v\f$ in the above equation refers to the valence bands of a
    semiconductor or an insulator.

    For detailed information see
    <a href="http://prola.aps.org/abstract/PRB/v59/i15/p10031_1"> Phys. Rev.
    B, 10031 (1999) </a>.

    \author Matthias Wahn, wahn@fhi-berlin.mpg.de
 */
class SX_EXPORT_EXX SxFockGk
{
   public:

      // --- all other than FCC and Wurtzite are not yet tested/reliable!!!
      enum Lattice { FCC, Wurtzite, SC, BCC, Cos, Gauss, None };
      
      /**@name Constructors and Destructor */
      SxFockGk ();
      //SxFockGk (const SxGkBasis &gk_, int nG_);
      SxFockGk (const SxGkBasis &gk_, int nG_, const SxCell &cell);
      virtual ~SxFockGk ()  { /* empty */ }

      void compute (int ik, const SxPW &waves, const SxFermi &fermi);

      void computeSymFock (int ik, const SxPW &waves, const SxFermi &fermi);

      PsiG operator* (const PsiG &psiG);

   protected:

      SxDiracMat<TPrecCoeffG>    Fock;

      SxDiracSymMat<TPrecCoeffG> symFock;

      Lattice                    lattice;

      SxCell                     cell;

      SxCell::CellType           cellType;

      const int                  nG;

      const SxGkBasis           *gkPtr;

      const int                  nk;

      const int                  nGkMax;

      SxVector<TPrecG>           kqG;

      void init ();

      void setupAuxTables ();

      void setupSingTable ();

      SxVector<TPrecG>           singTable;

      inline PrecG scF (const SxVector3<TPrecG> &p, PrecG a);
      inline PrecG fccF (const SxVector3<TPrecG> &p, PrecG a);
      inline PrecG bccF (const SxVector3<TPrecG> &p, PrecG a);
      inline PrecG wurtziteF (const SxVector3<TPrecG> &p, PrecG a, PrecG c);
      inline PrecG cosF (const SxVector3<TPrecG> &p, double pMax);
      inline PrecG gaussF (const SxVector3<TPrecG> &p, double lambda);

      /** Create table for the integral of the singularity function for
          the Wurtzite structure.
       */
      void createWurtziteTable ();

      /** sample points of the help function integral in case of
          the Wurtzite structure */
      SxDiracMat<TPrecG> wurtziteTable;

      /** Get the integral of the singularity help function over
          the Briouin zone by means of a linear interpolation between
          sampling points given in Ref. 2.

          return value of the integral in units of \f$(a/2\pi)^{2}/f$
       */
      PrecG interpolateWurtziteF (const PrecG caRatio);
};

#endif /* _SX_FOCK_GK_H_ */
