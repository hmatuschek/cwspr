#include "keyer.hh"

#define Nrise 64
#define dA 1./Nrise;


Keyer::Keyer(Source *source, QObject *parent)
  : Source(parent), _src(source), _state(STATE_OFF), _amp(0)
{
  if (_src) _src->setParent(this);
}

Keyer::~Keyer() {
  // pass...
}

void
Keyer::key(bool down) {
  if (down && (STATE_OFF == _state)) {
    _state = STATE_RISE;
    _amp = 0;
  } else if ((! down) && (STATE_ON == _state)) {
    _state = STATE_FALL;
    _amp = Nrise;
  }
}

void
Keyer::ditKey(bool down) {
  // pass...
}

void
Keyer::daKey(bool down) {
  // pass...
}

void
Keyer::read(int16_t *samples, qint64 nsamples) {
  if (_src)
    _src->read(samples, nsamples);
  else
    memset(samples, 0, sizeof(int16_t)*size_t(nsamples));

  for (int i=0; i<nsamples; i++) {
    samples[i] = int16_t((int32_t(samples[i])*_amp)/Nrise);
    if (STATE_RISE == _state) {
      if (Nrise == ++_amp)
        _state = STATE_ON;
    } else if (STATE_FALL == _state) {
      if (0 == --_amp)
        _state = STATE_OFF;
    }
  }

}
