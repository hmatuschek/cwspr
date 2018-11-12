#include "wsprencoder.hh"
#include <cmath>

WSPREncoder::WSPREncoder(float F0, const std::string &message)
    : Source(), _buffer(nullptr)
{
  setup(F0, message);
}

WSPREncoder::~WSPREncoder()
{
  if (_buffer)
    delete[] _buffer;
}

void
WSPREncoder::setup(float F0, const std::string &message, bool sync, int64_t phase) {
  _F0 = F0;
  _delay = ModeWSPR.delay*Fs;
  _phase = phase;
  _period = int(Fs)*int(ModeWSPR.period);
  _buffer = new int16_t[_period];
  int64_t symlen = int(Fs/ModeWSPR.sym_rate);
	std::vector<int> symbols;
	double phi=0;
	float wdt[4];

  WSPR::gen_wspr(message, symbols);
  for (int i=0; i<4; i++) {
    wdt[i] = 2*M_PI*(_F0+i*ModeWSPR.fsk_shift)/Fs;
  }

  for (int64_t phase=0; phase<_period; phase++) {
    _buffer[phase] = int16_t( (1<<15)*(std::sin(phi)) );
    size_t sidx = size_t(phase/symlen);
    if (sidx >= symbols.size())
      sidx = 0;
    phi += wdt[symbols[sidx]];
    if (phi >= 2*M_PI)
      phi -= 2*M_PI;
  }

  _sync = sync;
}

void
WSPREncoder::sync() {
  if (_sync)
    return;

  time_t clk;
  time(&clk);
  struct tm *utc = gmtime(&clk);
  int t = (60*utc->tm_min+utc->tm_sec) % int(_period/Fs);
  _phase = int64_t(t*Fs-_delay);
  if (_phase < 0)
    _phase = _period + _phase;
  _sync = true;
}


void
WSPREncoder::read(int16_t *samples, int64_t nsamples)
{
  sync();

  for (int64_t i=0; i<nsamples; i++) {
    samples[i] = _buffer[_phase++];
    if (_phase == _period)
      _phase = _phase - _period;
  }
}
