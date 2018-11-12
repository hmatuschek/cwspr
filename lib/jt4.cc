#include "jt4.hh"
#include <iostream>
#include <cstring>

Mode ModeJT4 = { MODE_JT4, 206, 60.0, 1.0, 4.375, 4.375 };


bool
encode4(const std::string &message, std::vector<int> &bits) {
  if (message.size() < 22)
    return false;
  bits.resize(207);
  encode4_(message.c_str(), bits.data());
  return true;
}


int
gen_jt4(const std::string &message, std::vector<int> &symbols) {
  int ichk=0, itype=0;
  char outm[22], msg[22];
  memset(msg, ' ', 22);
  memcpy(msg, message.c_str(), std::min(size_t(22), message.size()));
  symbols.resize(207);
  gen4_(msg, &ichk, outm, symbols.data(), &itype);
  if (ichk)
    return -ichk;
  return itype;
}
