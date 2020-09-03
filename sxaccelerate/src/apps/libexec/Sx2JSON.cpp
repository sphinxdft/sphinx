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

#include <SxCLI.h>
#include <SxSymbolTable.h>
#include <SxParser.tab.hpp>
#include <SxParser.h>
#include <SxFileIO.h>


void convertSymbol (const SxSymbol &in, SxList<SxString> *res,
                    const SxString &newLine, int lvl,
                    int indent, bool isStd)
{

   int t = in.type;
   switch (t)  {
      case VAR:  {
         // --- indent
         for (int j = 0; j < lvl*indent; ++j)  res->append (" ");

         if (isStd)  {
            if (in.name != "")  {
               if (in.name == "skipUnexpectedGroups")
                  res->append ("\"allowUndefined\":");
               else
                  res->append ("\"" + in.name + "\":");
            }

            if (   in.name == "optional"
                || in.name == "skipUnexpectedGroups"
                || in.name == "topLevelDefs")  {
               res->append ("true");
            }  else  {
               res->append (SxString(in.val));
            }
         }  else  {
            if (in.name != "")  res->append ("\"" + in.name + "\":");
            res->append (SxString(in.val));
         }
         break;
      }
      case STR:  {
         // --- indent
         for (int j = 0; j < lvl*indent; ++j)  res->append (" ");

         if (in.name != "")  res->append ("\"" + in.name + "\":");

         if (isStd && in.name == "type")  {
            if (in.str == "flag")       res->append ("\"bool\"");
            else                        res->append ("\"" + in.str + "\"");
         }  else  {
            res->append ("\"" + in.str + "\"");
         }

         break;
      }
      case LIST:  {
         // --- indent
         for (int j = 0; j < lvl*indent; ++j)  res->append (" ");

         if (in.name != "")  res->append ("\"" + in.name + "\":");
         ssize_t i = 0, n = in.valList->getSize ();
         res->append ("[");
         for (; i < n; ++i)  {
            res->append (newLine);
            convertSymbol ((*in.valList)(i), res, newLine, lvl+1, indent, isStd);
            if (i < n-1)  res->append (",");
         }
         res->append (newLine);
         // --- indent
         for (int j = 0; j < lvl*indent; ++j)  res->append (" ");

         res->append ("]");
         break;
      }
      default:
         SX_THROW ("invalid symbol type");
   }
}

void convertGroup (const SxSymbolTable *t, SxList<SxString> *res,
                   const SxString &newLine, int lvl,
                   int indent, bool isStd)
{
   SX_CHECK (t);
   SX_CHECK (res);

   // --- indent
   for (int j = 0; j < lvl*indent; ++j)  res->append (" ");

   if (t->name != "")  res->append ("\"" + t->getName () + "\":");

   res->append ("{");

   // --- convert all primitive type symbols
   SxMap<SxString, SxSymbol>::ConstIterator it = t->table.begin ();
   ssize_t i = 0, n = t->table.getKeys ().getSize ();

   for (; it != t->table.end (); ++it)  {
      res->append (newLine);
      convertSymbol (it.getValue (), res, newLine, lvl+1, indent, isStd);
      if (i < n-1)  res->append (",");
      ++i;
   }

   // --- convert all group type symbols

   // --- has printed any symbols
   bool hasSymbols = false;
   if (n > 0) hasSymbols = true;
   SxList<SxSymbolTable *>::ConstIterator lvlIt = t->children.begin ();
   i = 0;
   n = t->children.getSize ();

   // --- has printed symbols and has new children
   //     then add comma in between
   if (hasSymbols && n > 0)  res->append (",");

   for (; lvlIt != t->children.end (); ++lvlIt)  {
      res->append (newLine);
      convertGroup ((*lvlIt), res, newLine, lvl+1, indent, isStd);
      if (i < n-1)  res->append (",");
      ++i;
   }

   res->append (newLine);

   // --- indent
   for (int j = 0; j < lvl*indent; ++j)  res->append (" ");

   res->append ("}");
}


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   SxString inputFile  = cli.option ("--in", "file",
                                     "SPHIngX input file").toString ();

   SxString outputFile = cli.option ("--out", "file",
                                     "JSON output file").toString ();

   int indent = cli.option ("--tab", "indent",
                            "indentation").toInt (3);

   bool isCRNL = cli.option ("--crnl", "use windows style new line",
                             "use windows style new line").toBool ();

   bool isStd  = cli.option ("--std", "is std file",
                             "is std file").toBool ();

   SxString newLine = "\n";
   if (isCRNL == true)
      newLine = "\r\n";

   cli.finalize ();

   try  {

      SxParser parser;
      parser.setValidation (false);
      SxList<SxString> res;
      SxParser::Table table = parser.read (inputFile);
      res.append ("{");

      bool isFirst = true;
      if (isStd)  {
         // --- indent
         res.append (newLine);
         for (int j = 0; j < indent; ++j)  res.append (" ");
         res.append ("\"type\":\"group\"");
         isFirst = false;
      }
      // --- convert local variables
      for (auto tblIt = table->table.begin (); tblIt.isValid (); ++tblIt)  {
         if (   tblIt.getValue ().type == VAR
             || tblIt.getValue ().type == LIST
             || tblIt.getValue ().type == STR)
         {
            SxString name = tblIt.getValue ().name;
            bool invalid = (   name == "" || name == "PI"
                            || name == "Pi" || name == "pi"
                            || name == "true" || name == "TRUE"
                            || name == "false" || name == "FALSE"
                            || name == "yes" || name == "YES"
                            || name == "no" || name == "NO"
                            || name == "validator");

            if (!invalid)  {
               if (!isFirst)
                  res.append (",");
               res.append (newLine);
               convertSymbol (tblIt.getValue (), &res, newLine, 1, indent, isStd);
               isFirst = false;
            }
         }
      }

      // --- convert children
      SxList<SxSymbolTable *>::ConstIterator lvlIt = table->children.begin ();
      ssize_t i = 0;
      ssize_t n = table->children.getSize ();

      if (!isFirst && n > 0)  res.append (",");
      for (; lvlIt != table->children.end (); ++lvlIt)  {
         res.append (newLine);
         convertGroup ((*lvlIt), &res, newLine, 1, indent, isStd);
         if (i < n-1)  res.append (",");
         ++i;
      }

      res.append (newLine);

      res.append ("}");

      SxString resStr = SxString::join (res);

      SxFileIO::write (resStr, outputFile, 0600);

   }  catch (SxException e)  {
      e.print ();
   }
}
