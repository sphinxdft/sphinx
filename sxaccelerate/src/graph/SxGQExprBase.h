// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_GQ_EXPR_BASE_H_
#define _SX_GQ_EXPR_BASE_H_

#include <SxUniqueList.h>
#include <SxPtr.h>
#include <SxGProps.h>
#include <SxGraph.h>
#include <SxGQPattern.h>

/** \brief Graph Query Base Expression Class

    \b SxGQExprBase = SPHInX Graph Query Base Expression Class

     SxGQExprBase class represents the
     base class for all types of expressions.
 */
class SX_EXPORT_GRAPH  SxGQExprBase
{
   public:
      enum ExprType { None, Equal, NotEqual, Relation, Sibling, AND, OR, Any };

      ExprType exprType;

      typedef SxPtr<SxUniqueList<ssize_t> > Selection;
      typedef SxPtr<SxList<Selection> >  SelSet;

      SxGQExprBase() { }
      SxGQExprBase(ExprType t) : exprType(t) { }
      virtual ~SxGQExprBase() { }

      virtual SxPtr<SxList<SxPtr<SxGQExprBase> > > firsts () const = 0;
      virtual SxPtr<SxList<SxPtr<SxGQExprBase> > > lasts  () const = 0;
      virtual SxPtr<SxGQExprBase> first ()  const = 0;
      virtual SxPtr<SxGQExprBase> last  ()  const = 0;

      virtual void setLast  (const SxPtr<SxGQExprBase> &) = 0;
      virtual void setLasts (const SxPtr<SxList<SxPtr<SxGQExprBase> > > &) = 0;

      virtual bool       isOp ()       const = 0;
      virtual ExprType   getRightOp () const = 0;

      virtual size_t getHash () const = 0;

      virtual void makeGraph (SxGraph<SxGQPattern> *gPtr) const = 0;
      virtual SxGQPattern getGraphNode () const = 0;

      friend std::ostream &operator<< (std::ostream &s, const SxGQExprBase &item)
      {
         return item.print (s);
      }

   protected:
      virtual std::ostream &print (std::ostream &s) const = 0;
};
#endif /*_SX_GQ_EXPR_BASE_H_*/
