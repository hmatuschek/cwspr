#ifndef WSPRENCODER_HH
#define WSPRENCODER_HH

#include "wspr.hh"
#include "interfaces.hh"

class WSPREncoder: public WSPR, public Source
{
public:
  WSPREncoder(float F0, const std::string &message);
  virtual ~WSPREncoder();

  void setup(float F0, const std::string &message, bool sync=false, int64_t phase=0);

  void read(int16_t *samples, int64_t nsamples);

protected:
	void sync();

protected:
  float _F0;

  bool _sync;
  int64_t _delay;
  int64_t _phase;
  int64_t _period;
  int16_t *_buffer;

};

#endif // WSPRENCODER_HH
