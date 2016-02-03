#ifndef GOODRUN_H
#define GOODRUN_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <set>
#include <string>

struct run_and_lumi{
  unsigned int run;
  long long int lumi_min;
  long long int lumi_max;
};

bool operator<(const struct run_and_lumi &r1, const struct run_and_lumi &r2)
{
  return r1.run < r2.run;
}


class GoodRun {
  
  enum file_type{TEXT, JSON};

  typedef std::multiset<struct run_and_lumi> set_t; 
  set_t good_runs_;
  bool good_runs_loaded_ = false;
  
 public:
  int load_runs(const char *fname, enum file_type type);
  bool goodrun (unsigned int run, unsigned int lumi_block);
  bool goodrun_json (unsigned int run, unsigned int lumi_block);
  void set_goodrun_file (const char* filename);
  void set_goodrun_file_json (const char* filename);
  int  min_run ();
  int  max_run ();
  int  min_run_min_lumi ();
  int  max_run_max_lumi ();

};

#endif
