#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <string.h>
#include <errno.h>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include "h_bond.h"
#include "Fresidue.h"
#include "control.h"
#include "element.h"
#include "mpi.h"
using namespace std;

// Command line. Everything here is argument handling; the matrix itself is
// built by control.h exactly as before, so records are byte-identical to the
// v1.0.0 build (checked by scripts/db_validation/genmat_args_regression.sh).
//
// The original front end failed silently in ways a reasonable caller hits:
// a directory without a trailing slash was concatenated straight onto the file
// name, a missing input exited 0 after creating an empty output file, and the
// batch modes -ds/-fs hang without mpirun and write zero-element records under
// it (the worker path never calls prepareIndex). All of those are now errors.

static const char *USAGE =
"usage: generateMatrix [-os|-o] <file.ssd> <output_file>\n"
"       generateMatrix [-os|-o] <name.ssd> <palsse_dir> <output_file>\n"
"\n"
"Builds the interaction matrix for ONE PALSSE .ssd file and writes it to\n"
"<output_file>. Concatenate the outputs into a searchable DB:\n"
"  cat out/*.out > metamatricesDB\n"
"\n"
"  -os  sheets from the PALSSE file's SHEET records (default; what the\n"
"       census DBs were built with)\n"
"  -o   sheets recomputed by generateMatrix (legacy)\n"
"\n"
"The record is named after the .ssd file's basename. No mpirun is needed:\n"
"run one process per file and parallelise with xargs -P or a SLURM array.\n"
"The old batch modes -ds and -fs are disabled; they hang, or write records\n"
"with zero elements.\n";

static int finish(int code)
{
  MPI_Finalize();
  return code;
}

static string joinPath(const string &dir, const string &name)
{
  if(dir.empty()) return name;
  if(dir[dir.size()-1] == '/') return dir + name;
  return dir + "/" + name;
}

static string baseName(const string &path)
{
  size_t slash = path.find_last_of('/');
  return slash == string::npos ? path : path.substr(slash + 1);
}

// control's readers take char*; give them a writable, unbounded copy.
static vector<char> cbuf(const string &s)
{
  vector<char> b(s.begin(), s.end());
  b.push_back('\0');
  return b;
}

int main(int argc, char* argv[])
{
  int numprocs, myrank;
  control con;
  vector<element> eleInOneStru;

  MPI_Init(&argc ,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank != 0)
     return finish(0);   // single-structure tool: extra ranks have nothing to do

  vector<string> args(argv + 1, argv + argc);
  if(args.empty())
  {
     fputs(USAGE, stderr);
     return finish(2);
  }
  if(args[0] == "-h" || args[0] == "--help")
  {
     fputs(USAGE, stdout);
     return finish(0);
  }

  string option = "-os";
  if(args[0][0] == '-')
  {
     option = args[0];
     args.erase(args.begin());
     if(option == "-ds" || option == "-d" || option == "-fs" || option == "-f")
     {
        fprintf(stderr, "generateMatrix: %s is disabled: it hangs without mpirun and "
                        "writes zero-element records under it.\n"
                        "Run generateMatrix -os once per .ssd file instead.\n\n", option.c_str());
        fputs(USAGE, stderr);
        return finish(2);
     }
     if(option != "-os" && option != "-o")
     {
        fprintf(stderr, "generateMatrix: unknown option %s\n\n", option.c_str());
        fputs(USAGE, stderr);
        return finish(2);
     }
  }

  string ssdPath, name, outPath;
  if(args.size() == 2)
  {
     ssdPath = args[0];
     outPath = args[1];
  }
  else if(args.size() == 3)
  {
     ssdPath = joinPath(args[1], args[0]);
     outPath = args[2];
  }
  else
  {
     fprintf(stderr, "generateMatrix: expected 2 or 3 arguments after the option, got %d\n",
             (int)args.size());
     for(size_t i = 1; i < args.size(); i++)
        if(args[i][0] == '-')
        {
           fprintf(stderr, "(%s must come first: generateMatrix %s <file.ssd> <output_file>)\n",
                   args[i].c_str(), args[i].c_str());
           break;
        }
     fputs("\n", stderr);
     fputs(USAGE, stderr);
     return finish(2);
  }
  name = baseName(ssdPath);

  // Check the input BEFORE creating the output, so a bad path leaves no empty
  // file behind to be concatenated into a DB.
  FILE *probe = fopen(ssdPath.c_str(), "r");
  if(probe == NULL)
  {
     fprintf(stderr, "generateMatrix: cannot read %s: %s\n", ssdPath.c_str(), strerror(errno));
     return finish(1);
  }
  fclose(probe);

  FILE *a = fopen(outPath.c_str(), "w");
  if(a == NULL)
  {
     fprintf(stderr, "generateMatrix: cannot write %s: %s\n", outPath.c_str(), strerror(errno));
     return finish(1);
  }

  vector<char> fileName = cbuf(ssdPath);
  vector<char> tokenstring = cbuf(name);
  con.setoption((char*)option.c_str());
  cout<<"the file name is "<<&fileName[0]<<endl;
  con.prepareIndex(&fileName[0]);
  con.readIndrefile1(&fileName[0], eleInOneStru);
  cout<<"this out of read file"<<endl;
  con.Elmentdelete(eleInOneStru);
  con.prodvector(eleInOneStru);
  if(eleInOneStru.empty())
     fprintf(stderr, "generateMatrix: warning: %s has no helix or strand elements; "
                     "writing an empty record\n", ssdPath.c_str());
  con.proInteractionMatr(eleInOneStru, a, &tokenstring[0]);
  if(fclose(a) != 0)
  {
     fprintf(stderr, "generateMatrix: error writing %s: %s\n", outPath.c_str(), strerror(errno));
     return finish(1);
  }
  return finish(0);
}
