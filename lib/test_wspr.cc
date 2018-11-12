#include <iostream>

#include <cinttypes>
#include "cwspr.hh"
#include <unistd.h>


int main(int argc, char *argv[])
{
  int64_t N = 12000*120;
  int16_t *buffer = new int16_t[N];

  WSPREncoder encoder(600., "DM3MAT JO62 37");
  WSPRDecoder decoder(nullptr, 600.);

  encoder.setup(600., "DM3MAT JO62 37", true, 0);
  decoder.setup(600, true, 0);

  encoder.read(buffer, N);
  decoder.write(buffer, N);

  sleep(2);
  return 0;
}
