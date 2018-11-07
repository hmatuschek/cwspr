#ifndef CWSPR_WSPR_HH
#define CWSPR_WSPR_HH

#include <string>
#include <vector>
#include "traits.hh"


extern "C" {
 void wqencode_(const char *msg, int *ntype, char *data);
}

int gen_wspr(const std::string &message, std::vector<int> &symbols);

extern Mode ModeWSPR;

#endif // CWSPR_WSPR_HH
