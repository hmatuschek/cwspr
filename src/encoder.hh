#ifndef ENCODER_HH
#define ENCODER_HH

#include "source.hh"
#include "traits.hh"

class Encoder : public Source
{
	Q_OBJECT

public:
	Encoder(ModeId modeId, const QString &message, double freq, QObject *parent=nullptr);
	virtual ~Encoder();

	bool setup(ModeId modeId, const QString &message, double freq);

	virtual void read(int16_t *samples, qint64 nsamples);

protected:
	void sync();

protected:
	bool _sync;
	qint64 _delay;
	qint64 _phase;
	qint64 _period;
	int16_t *_buffer;
};


#endif // ENCODER_HH
