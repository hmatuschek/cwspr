#ifndef CWSPR_APPLICATION_HH
#define CWSPR_APPLICATION_HH

#include <QApplication>
#include <portaudio.h>
#include "encoder.hh"
#include "keyer.hh"
#include "waterfall.hh"


class Application : public QApplication
{
	Q_OBJECT

public:
	Application(int &argc, char *argv[]);
	virtual ~Application();

	void key(bool down);
	void ditKey(bool down);
	void daKey(bool down);

	WaterfallSink *waterfall();

public slots:
	void start();
	void stop();

protected slots:
	void info();

protected:
	static int _pa_outstream_callback(const void *in, void *out, unsigned long frameCount,
	                                  const PaStreamCallbackTimeInfo *tInfo,
	                                  PaStreamCallbackFlags status, void *userData);
	static int _pa_instream_callback(const void *in, void *out, unsigned long frameCount,
	                                 const PaStreamCallbackTimeInfo *tInfo,
	                                 PaStreamCallbackFlags status, void *userData);

protected:
	PaStream *_out;
	PaStream *_in;
	Encoder *_encoder;
	Keyer *_keyer;
	WaterfallSink *_waterfall;
};

#endif // CWSPR_APPLICATION_HH
