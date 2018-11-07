#include "jt4.hh"

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
gen_jt4(const std::string &message, std::vector<int> &symbols, std::string &outmsg) {
  int ichk, itype;
  char outm[22];
  symbols.resize(207);
  gen4_(message.c_str(), &ichk, outm, symbols.data(), &itype);
  outmsg = outm;
  if (ichk)
    return -ichk;
  return itype;
}
