#ifndef KEYER_HH
#define KEYER_HH

#include "source.hh"


class Keyer: public Source
{
	Q_OBJECT

private:
	typedef enum { STATE_OFF, STATE_RISE, STATE_ON, STATE_FALL } State;

public:
  Keyer(Source *source, QObject *parent=nullptr);
  virtual ~Keyer();

  Source *source() const;
  void setSource(Source *source);

  void read(int16_t *samples, qint64 nsamples);

  void key(bool down);
  void ditKey(bool down);
  void daKey(bool down);

protected:
  Source *_src;
  State _state;
  int16_t _amp;
};

#endif // KEYER_HH
