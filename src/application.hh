#ifndef CWSPR_APPLICATION_HH
#define CWSPR_APPLICATION_HH

#include <QApplication>
#include <portaudio.h>
#include "encoder.hh"
#include "keyer.hh"
#include "waterfall.hh"
#include <QTimer>
#include "wsprdecoder.hh"


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

signals:
	void ptt(bool tx);

protected slots:
	void info();
	void onPTTTimeOut();

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
  WSPRDecoder *_decoder;

	Keyer *_keyer;
	WaterfallSink *_waterfall;
	bool _ptt;
	QTimer _pttTimer;
};

#endif // CWSPR_APPLICATION_HH
