#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#ifdef SILENT
#include <cstdio>
#endif
#include "elecol.h"
#include "searchControl.h"
#include "Fpass.h"
#include "external.h"
#include "MtrixElment.h"
#include "handness.h"
#include "sheet.h"
#include "require.h"
using namespace std;

// Command line. Everything here is argument checking; the search itself is
// searchControl::oneprocess, unchanged, so hit output is byte-identical to the
// v1.0.0 build (scripts/db_validation/searchmatrix_args_regression.sh).
//
// With -DSILENT (the Makefile build) stdout is /dev/null, so the old usage
// text and every "can't open" message vanished and the run exited 0 with no
// hits: a missing DB, a missing query, a manifest with a bad path, too few
// arguments. Those are now checked here, before stdout is silenced, and
// reported on stderr with a nonzero exit.

static const char *USAGE =
"usage: searchmatrix <query.query | manifest> <metamatricesDB> <output_dir>\n"
"\n"
"  query     one ProSMoS query file; hits go straight into <output_dir>\n"
"  manifest  a file listing query files, one path per line; hits for each go\n"
"            into <output_dir>/<query name>/ (the DB is read once for all)\n"
"  output_dir  created if missing; the trailing slash is optional\n"
"\n"
"One hit file per matching DB record. Search the .clean DB (see\n"
"scripts/db_validate.py). Exit status: 0 ok, 1 unreadable input or a\n"
"malformed query/DB, 2 bad arguments. Query format: README, \"Query format\".\n";

static bool readable_file(const string &p)
{
  struct stat st;
  return stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode) && access(p.c_str(), R_OK) == 0;
}

static bool make_dirs(const string &path)
{
  string cur;
  size_t pos = 0;
  while(pos != string::npos)
  {
    size_t next = path.find('/', pos + 1);
    cur = path.substr(0, next);
    if(!cur.empty() && mkdir(cur.c_str(), 0777) != 0 && errno != EEXIST)
      return false;
    pos = next;
  }
  struct stat st;
  return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

static bool query_number_line(const string &l)
{
  if(l.empty()) return false;
  for(size_t i = 0; i < l.size(); i++)
    if(!isdigit((unsigned char)l[i]) && l[i] != ' ' && l[i] != '\t' && l[i] != '\r')
      return false;
  return true;
}

int main(int argc , char *argv[])
{
  vector<string> args(argv + 1, argv + argc);
  if(args.size() == 1 && (args[0] == "-h" || args[0] == "--help"))
  {
     fputs(USAGE, stdout);
     return 0;
  }
  if(args.size() != 3)
  {
     fprintf(stderr, "searchmatrix: expected 3 arguments, got %d\n\n", (int)args.size());
     fputs(USAGE, stderr);
     return 2;
  }
  string qarg = args[0], db = args[1], out = args[2];

  // A query file starts with its element-number line ("1 2 3 4"); anything
  // else is a manifest, and then every listed query must exist. (The engine
  // guesses by stat()ing the first line, so a manifest whose first path was
  // wrong used to be parsed as a query.)
  if(!readable_file(qarg))
  {
     fprintf(stderr, "searchmatrix: cannot read %s: %s\n", qarg.c_str(),
             access(qarg.c_str(), F_OK) == 0 ? "not a readable file" : strerror(errno));
     return 1;
  }
  {
     ifstream in(qarg.c_str());
     string first;
     getline(in, first);
     if(!query_number_line(first))
     {
        ifstream mf(qarg.c_str());
        string line;
        int n = 0, bad = 0, lineno = 0;
        while(getline(mf, line))
        {
           lineno++;
           while(!line.empty() && isspace((unsigned char)line[line.size()-1])) line.erase(line.size()-1);
           if(line.empty())
           {
              if(lineno == 1)
              {
                 fprintf(stderr, "searchmatrix: %s: first line is empty (neither a query nor a manifest)\n", qarg.c_str());
                 return 1;
              }
              continue;
           }
           n++;
           if(!readable_file(line))
           {
              if(bad < 5)
                 fprintf(stderr, "searchmatrix: manifest %s line %d: cannot read query '%s'\n",
                         qarg.c_str(), lineno, line.c_str());
              bad++;
           }
        }
        if(bad)
        {
           fprintf(stderr, "searchmatrix: %d of %d manifest entries unreadable%s\n", bad, n,
                   bad > 5 ? " (first 5 shown)" : "");
           return 1;
        }
        if(n == 0)
        {
           fprintf(stderr, "searchmatrix: manifest %s lists no queries\n", qarg.c_str());
           return 1;
        }
     }
  }

  if(!readable_file(db))
  {
     fprintf(stderr, "searchmatrix: cannot read DB %s: %s\n", db.c_str(),
             access(db.c_str(), F_OK) == 0 ? "not a readable file" : strerror(errno));
     return 1;
  }

  while(out.size() > 1 && out[out.size()-1] == '/') out.erase(out.size()-1);
  if(!make_dirs(out) || access(out.c_str(), W_OK) != 0)
  {
     fprintf(stderr, "searchmatrix: cannot create or write output dir %s: %s\n", out.c_str(), strerror(errno));
     return 1;
  }
  out += "/";

#ifdef SILENT
  // searchControl.h has 200+ leftover debug cout/printf calls inside the
  // per-DB-entry match loop; on a 710k-entry DB that's tens of millions of
  // small writes per query. Redirect stdout once at startup -- both cout and
  // printf hit fd 1, so this catches them all without touching the call sites.
  // Errors go to cerr (fd 2) and exit nonzero, so they still surface.
  FILE *silent_fp = freopen("/dev/null", "w", stdout);
  (void)silent_fp;
#endif
  searchControl control;
  vector<char> a1(qarg.begin(), qarg.end()); a1.push_back('\0');
  vector<char> a2(db.begin(), db.end());     a2.push_back('\0');
  vector<char> a3(out.begin(), out.end());   a3.push_back('\0');
  cout<<"this is step one : "<<endl;
  control.oneprocess(&a1[0], &a2[0], &a3[0]);
  return 0;
}
