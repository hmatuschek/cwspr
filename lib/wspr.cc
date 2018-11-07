#include "wspr.hh"
#include <iostream>


extern "C" {
void encode232_(const char *data, int *dlen, char *sym);
void inter_wspr_(char *sym, int *dir);
}

const int sync_[162] = {
  1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
  0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
  0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
  1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
  0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
  0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
  0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
  0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
  0,0 };

Mode ModeWSPR = { MODE_WSPR, 162, 120.0, 1.0, 1.464738, 1.4648};

int
gen_wspr(const std::string &message, std::vector<int> &symbols) {
  char msg[22], data[13];
  memset(msg, ' ', 22);
  memcpy(msg, message.c_str(), std::min(size_t(21), message.size()));
  memset(data, 0, 13);

  int ntype=0;
  wqencode_(msg, &ntype, data);
  std::cerr << "Type " << ntype << std::endl;
  for (int i=0; i<13; i++)
    std::cerr << " " << std::hex << int(data[i]);
  std::cerr << std::endl;

  char sym[206];
  memset(sym, 0, 206);
  int slen=206, dir=1;
  encode232_(data, &slen, sym);
  inter_wspr_(sym, &dir);

  symbols.resize(162);
  std::cerr << "Sym:";
  for (int i=0; i<162; i++) {
    symbols[i] = 2*int(sym[i])+sync_[i];
    std::cerr << " " << int(symbols[i]);
  }

  return ntype;
}
