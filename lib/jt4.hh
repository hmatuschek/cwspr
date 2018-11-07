#ifndef CWSPR_JT4_H
#define CWSPR_JT4_H

#include "traits.hh"
#include <string>
#include <vector>

extern "C" {
// interface to Fortran impl.
void encode4_(const char *msg, int *sym);
// interface to Fortran impl.
void gen4_(const char *msg, int *ichk, char *msgsent, int *sym, int *itype);
}


/** Encodes the given message.
 */
bool encode4(const std::string &message, std::vector<int> &bits);

/** Generates symbol codes for the given message.
 */
int gen_jt4(const std::string &message, std::vector<int> &symbols);


extern Mode ModeJT4;

#endif // CWSPR_JT4_H
