#include "encoder.hh"
#include "jt4.hh"
#include "wspr.hh"
#include <QDebug>
#include <cmath>

#define Fs 16000.

Encoder::Encoder(ModeId mode, const QString &message, double freq, QObject *parent)
  : Source(parent), _sync(false), _phase(0), _period(0), _buffer(nullptr)
{
  setup(mode, message, freq);
}

Encoder::~Encoder() {
  if (_buffer)
    delete [] _buffer;
}

bool
Encoder::setup(ModeId modeId, const QString &message, double freq)
{
  Mode mode = (MODE_WSPR == modeId) ? ModeWSPR : ModeJT4;
  qDebug() << "setup mode " << modeId;

  _delay = mode.delay*Fs;
  _phase = 0;
  _period = int(Fs)*int(mode.period);
  _buffer = new int16_t[_period];
  qint64 symlen = int(Fs/mode.sym_rate);
	std::vector<int> symbols;
	double phi=0;
	float wdt[4];

  std::string omsg;
  if (MODE_WSPR == modeId)
    gen_wspr(message.toStdString(), symbols);
  if (MODE_JT4 == modeId)
    gen_jt4(message.toStdString(), symbols, omsg);

  for (int i=0; i<4; i++) {
    wdt[i] = 2*M_PI*(freq+i*mode.fsk_shift)/Fs;
  }

  for (qint64 phase=0; phase<_period; phase++) {
    _buffer[phase] = int16_t( (1<<15)*(std::sin(phi)) );
    size_t sidx = size_t(phase/symlen);
    if (sidx >= symbols.size())
      sidx = 0;
    phi += wdt[symbols[sidx]];
    if (phi >= 2*M_PI)
      phi -= 2*M_PI;
  }

  _sync = false;

  return true;
}


void
Encoder::sync() {
  if (_sync)
    return;

  time_t clk;
  time(&clk);
  struct tm *utc = gmtime(&clk);
  int t = (60*utc->tm_min+utc->tm_sec) % int(_period/Fs);
  _phase = qint64(t*Fs-_delay);
  if (_phase < 0)
    _phase = _period + _phase;
  _sync = true;
}


void
Encoder::read(int16_t *samples, qint64 nsamples) {
  sync();

  for (qint64 i=0; i<nsamples; i++) {
    samples[i] = _buffer[_phase++];
    if (_phase == _period)
      _phase = _phase - _period;
  }
}

